# -*- coding: utf-8 -*-
""" CRG-QC2-View: `panel serve qc2_timeseries.py` """
import os
import datetime as dt
from typing import Dict, List
from dataclasses import dataclass, field
import logging
import json

import pandas as pd
from mzqc import MZQCFile as qc

import panel as pn
import hvplot.pandas
import holoviews as hv

pn.extension('tabulator')
hv.extension('bokeh')
pn.extension(loading_spinner='dots', loading_color='#633663', template='bootstrap')

qc2log = logging.getLogger()

DESCR_TEMPL = """
# The Data Set

This dataset ({p}) is comprised of mzQC files from {n} runs,

which were measured over the timecourse between {s} and {e}, 

from QC2 samples on an instrument of type {i}.

You can in inspect the details in the following interactive charts.
"""

@dataclass
class dataset:
  """Class for keeping track the different metric dataframes for a particular dataset.      
      Adds a "dRT" column to the peps dataframe which is the retention time delta of the 
      respective peptide in a given run from the datasets mean retention time of that peptide.
      Updates the dataset object with further statistical information on the QC peptides in two 
      dicts (for the mean and the standard deviation of the respective) values of:
        "dRT": deviations from the whole-dataset mean of the respective peptides' retention time,
        "dPPM": mass deviations of QC peptides found,
        "# Chargestates": number of unique charge states in which the QC peptides were found per run,
        "# Identified QC2 Peptides": number of QC2 peptides identified per run,
  """
  main: pd.DataFrame
  tics: pd.DataFrame 
  peps: pd.DataFrame
  meta: pd.DataFrame
  label: str
  instrument: str
  peps_std: Dict[str,float] = field(default_factory=dict)
  peps_mean: Dict[str,float] = field(default_factory=dict)

  def start_str(self):
    return str(self.main.Date.min().date())
  def end_str(self):
    return str(self.main.Date.max().date())
  def size_str(self):
    return str(len(self.main))
  def main_index_to_names(self, selection_indices):
    return self.main.Name.iloc[selection_indices].to_list()
  
  def __post_init__(self):
    self.peps = pd.merge(self.peps, pd.DataFrame(
                  {"mean_rt_pp": self.peps.groupby(['Peptide'])['RT'].mean()}), how="inner", on="Peptide")
    self.peps["dRT"] = (self.peps['mean_rt_pp'] - self.peps['RT']).drop(columns=["mean_rt_pp"])
    self.peps_std = {
      "dRT": self.peps['dRT'].std(),
      "dPPM": self.peps['dPPM'].std(),
      "# Chargestates": self.peps.groupby(['Date','Name'])['Chargestate'].nunique().std(),
      "# Identified QC2 Peptides": self.peps.groupby(['Date','Name'])['Peptide'].nunique().std(),
    }
    self.peps_mean = {
      "dRT": self.peps['dRT'].mean(),
      "dPPM": self.peps['dPPM'].mean(),
      "# Chargestates": self.peps.groupby(['Date','Name'])['Chargestate'].nunique().mean(),
      "# Identified QC2 Peptides": self.peps.groupby(['Date','Name'])['Peptide'].nunique().mean(),
    }

def extract_date(rowrunmzqc: qc.MzQcFile) -> pd.Timestamp:
  """extracts a time object from a mzQC object

  Parameters
  ----------
  rowrunmzqc : qc.MzQcFile
      the mzQC object to load the run date from (expects "MS:1000747" to be present)

  Returns
  -------
  pd.Timestamp
      a pandas Timestamp; takes the first run found in the mzQC object.

      Picks the date source data from the completion time from the first run in the mzQC object.
      TODO pick the run by provided name and fail gracefully if not found.

  """
  rawdate = next(iter([pd.to_datetime(cd.value) for cd in rowrunmzqc.runQualities[0].metadata.inputFiles[0].fileProperties if cd.accession=="MS:1000747"]))
  return pd.Timestamp(rawdate)

def load_main_from_mzqc(rowrunname: str, rowrunmzqc: qc.MzQcFile) -> pd.DataFrame:
  """loads a single values QC metric dataframe from a mzQC object

  Parameters
  ----------
  rowrunname : str
      the run name for the contents of the respective mzQC object
  rowrunmzqc : qc.MzQcFile
      the mzQC object to load the tics from (expects "MS:4000059","MS:4000060","MS:1003251",
      "MS:1003250","MS:1003328","MS:1002404","MS:4000102","MS:4000103","MS:4000097","MS:4000098",
      "MS:4000070","MS:4000069","MS:4000069" to be present)

      Takes the first run found in the mzQC object.
      TODO pick the run by provided name and fail gracefully if not found.

  Returns
  -------
  pd.DataFrame
      single row dataframe; columns=['# MS1','# MS2','# ID MS2','# Peptidoforms',
      '# Proteoforms','# Proteins','# Features','# ID Features',u'# Signal fluct. ↑',
      u'# Signal fluct. ↓','RT range left','RT range right','MZ range left',
      'MZ range right','Name'

  See Also
  --------
  extract_date : Defines the type of the Date column.
  """
  rundate = extract_date(rowrunmzqc)
  df = pd.DataFrame({'Date': rundate,
          'Name':rowrunname,
          'Date .. Name': rundate.strftime('%Y-%m-%d') + ' .. ' + rowrunname[-10:],
          '# MS1': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000059"])),
          '# MS2': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000060"])),
          '# ID MS2': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:1003251"])),
          '# Peptidoforms': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:1003250"])),
          '# Proteoforms': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:1003328"])),
          '# Proteins': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:1002404"])),
          '# Features': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000102"])),
          '# ID Features': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000103"])),
          u'# Signal fluct. ↑': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000097"])),
          u'# Signal fluct. ↓': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000098"])),
          # 'RT range': str(next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000070"]))),
          # 'MZ range': str(next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000069"]))),
          'mzrange': [cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000069"],
          'rtrange': [cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000070"],
          })
  df = df.join(pd.DataFrame(df.pop('mzrange').tolist(), index=df.index, columns=["MZ range left", "MZ range right"]))
  df = df.join(pd.DataFrame(df.pop('rtrange').tolist(), index=df.index, columns=["RT range left", "RT range right"]))

  return df

def load_tic_from_mzqc(rowrunname: str, rowrunmzqc: qc.MzQcFile) -> pd.DataFrame:
  """loads a tic dataframe from a mzQC file object

  Parameters
  ----------
  rowrunname : str
      the run name for the contents of the respective mzQC object
  rowrunmzqc : qc.MzQcFile
      the mzQC files object to load the tics from (expects "MS:4000104" to be present)

      Creates the obvious cols of "RT"("MS:1000894"), "Intensity" ("MS:1000285"), "# Peak"("MS:1003059"), "NativeID"("MS:1000767"),
      but also 3 convenience cols: "Date" (taken from the mzQC), "Name" (as provided by the fn param, to keep the data adressable 
      after merging with all the other tic dataframes, same for Date), and "Date .. Name" as contraction of the date and the last 
      10 characters of the run name for convenient legend display if there are more than one run on the date.
      
      Takes the first run found in the mzQC.
      TODO pick the run by provided name and fail gracefully if not found

  Returns
  -------
  pd.DataFrame 
      dataframe as produced above; columns=['RT', 'Intensity', 'Date', 'Name', 'Date .. Name']

  See Also
  --------
  extract_date : Defines the type of the Date column.
  """
  date = extract_date(rowrunmzqc)
  tics = pd.DataFrame(next(iter([m.value for m in rowrunmzqc.runQualities[0].qualityMetrics if m.accession=="MS:4000104"])))\
              .rename(columns={"MS:1000894": "RT", "MS:1000285": "Intensity", "MS:1003059": "# Peak", "MS:1000767": "NativeID"})\
              .assign(**{"Date":date, "Name":rowrunname, "Date .. Name": date.strftime('%Y-%m-%d') + " .. " + rowrunname[-10:]})
  return tics
  
def load_peps_from_mzqc(rowrunname: str, rowrunmzqc: qc.MzQcFile) -> pd.DataFrame:
  """loads a QC peptides dataframe from a mzQC file object

  Parameters
  ----------
  rowrunname : str
      the run name for the contents of the respective mzQC object
  rowrunmzqc : qc.MzQcFile
      the mzQC files object to load the tics from (expects "MS:4000078" to be present)
  
      Creates the obvious cols of "RT"("MS:1000894"), "Peptide" ("MS:1003169"), "dPPM"("MS:1003059"), "Chargestate"("MS:1000041"),
      but also 2 convenience cols: "Date" (taken from the mzQC), "Name" (as provided by the fn param, to keep the data adressable 
      after merging with all the other pep dataframes, same for Date).
      
      Takes the first run found in the mzQC.
      TODO pick the run by provided name and fail gracefully if not found

  Returns
  -------
  pd.DataFrame
      dataframe as produced above; columns=['RT','Peptide', 'dPPM', '# Chargestates', 'Date', 'Name']

  See Also
  --------
  extract_date : Defines the type of the Date column.
  """
  date = extract_date(rowrunmzqc)
  df = pd.DataFrame(next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000078"])))\
              .rename(columns={"MS:1000894": "RT",
                                "MS:1003169": "Peptide",
                                "MS:4000072": "dPPM",
                                "MS:1000041": "Chargestate"})\
              .assign(**{'Date': date, 'Name': rowrunname})
  df["Chargestate"] = df["Chargestate"].astype('int')
  df[["RT", "dPPM"]] = df[["RT", "dPPM"]].astype('float')
  return df

def load_ds(ds_key: str, mzqc_paths: Dict[str,List[str]], metadata_paths: Dict[str,List[str]]) -> dataset:
  """loads all the mzQC files from a given project and creates a dataset representation

  Parameters
  ----------
  ds_key : str
      the dataset key to the mzqc_paths parameter
  mzqc_paths : Dict[str,List[str]]
      the dict of lists for all mzQC file paths from the respective dataset
  
  Returns
  -------
  dataset
      the dataset object created from the mzQC files

      All dataframes are sorted by Date.
      
  See Also
  --------
  load_main_from_mzqc : Defines the dataframe for single value metrics.
  load_tic_from_mzqc : Defines the dataframe for ion chromatography metrics.
  load_peps_from_mzqc : Defines the dataframe for QC2 peptide metrics.
  """
  # adds 'dRT','# Identified QC2 Peptides' to pep?
  if ds_key not in mzqc_paths.keys():
    logging.error("ds_key not in mzqc_paths.")
    return
  qc2log.debug("updating with {}".format(ds_key))

  dfs_main = list()
  dfs_tics = list()
  dfs_peps = list()
  for mzqc_path in mzqc_paths.get(ds_key, None):
    with open(mzqc_path, "r") as file:
        mzqcobj = qc.JsonSerialisable.FromJson(file)
        name = os.path.basename(mzqc_path)
        dfs_main.append(load_main_from_mzqc(name,mzqcobj))
        dfs_tics.append(load_tic_from_mzqc(name, mzqcobj))
        dfs_peps.append(load_peps_from_mzqc(name,mzqcobj))

  with open(metadata_paths.get(ds_key, None), "r") as file:
    df_meta = pd.DataFrame(json.load(file)).\
      rename(columns={"user_date": "Date", "user_email": "Contact", 
                      "additional_information": "Additional_Information", 
                      "problems": "Problems",
                      "actions": "Actions"})
    df_meta.Date = df_meta.Date.astype('datetime64[ns]')
    df_meta.Actions = df_meta.Actions.apply(lambda x: '&'.join([y.get('name',None) for y in x]))

  return dataset(main = pd.concat(dfs_main, axis=0, ignore_index=True).sort_values(by='Date'), 
                peps = pd.concat(dfs_peps, axis=0, ignore_index=True).sort_values(by='Date'),
                tics = pd.concat(dfs_tics, axis=0, ignore_index=True).sort_values(by='Date'),
                meta =  df_meta.sort_values(by='Date'),
                label = dataset_mid_path[ds_key], 
                instrument = "Thermo {}".format(next(iter(ds_key.split('_')))))

def plot_calendar_hist(data:pd.DataFrame) -> hvplot.hvPlot:
  """plots a histogram of the number of runs, binned by month

  Parameters
  ----------
  data : pd.DataFrame
      dataframe, representing one run per row, containing at least a
      'Date' column for the run completion date (in datetime64['ns'])

  Returns
  -------
  hvplot.hvPlot
      a Holoviews plot object
  """
  qc2log.debug("plot_calendar histogram")
  if (data.shape[0] == 0):
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(
          title="QC2 Run Yearly Distribution", 
          ylabel="# Runs",
          xlabel="Month of {}".format("N/A"))

  qc2log.debug("plot_calendar histogram in ernest")
  calendar_hist = data.assign(month=data.Date.dt.month)\
        .groupby(by=['month']).count().hvplot.bar(y="Date")\
        .opts(title="QC2 Run Yearly Distribution", ylabel="# Runs",
              xlabel="Month of {}".format(next(iter(data.Date.dt.year)))).opts(shared_axes=False)
  return calendar_hist

def plot_tics(selection:List[str], data:pd.DataFrame) -> hvplot.hvPlot:
  """plots the tic from a 

  Parameters
  ----------
  selection : List[int]
      list of names of runs(, selected from the df_widget)
  data : pd.DataFrame
      the dataframe containing the tics of multiple runs in long format (columns=["Name","RT","Intensity","Date .. Name"])

  Returns
  -------
  hvplot.hvPlot
      a Holoviews plot object
  """
  qc2log.debug("plot_tics")
  if (data.shape[0]==0) or\
     (len(selection)==0):
    return pd.DataFrame(columns=["RT", "Intensity"]).hvplot.line(title="Total Ion Chromatogram",
                                                            line_width=0.5, 
                                                            xlabel="Retentiontime [s]", 
                                                            ylabel="Relative Intensity of Ion Current",
                                                            frame_height=500,
                                                            frame_width=1000,) * hv.Text(.5,.5,"Select Run\nfrom Table", fontsize=24,).opts(color='grey')
  
  qc2log.debug("plot_tics in ernest")
  # with pn.param.set_values(col[3][0], loading=True):
  while True:
    selected_tic = data[data.Name.isin(selection)].copy()
    selected_tic = selected_tic[["RT","Intensity","Date .. Name"]].sort_values(["RT"], axis = 0, ascending = True)
    plot = selected_tic.hvplot.line(x="RT",y="Intensity", by=["Date .. Name"],
                                    title="Total Ion Chromatogram",
                                    xlabel="Retentiontime [s]", ylabel="Relative Intensity of Ion Current",
                                    frame_height=500, frame_width=1000,).opts(shared_axes=False)
    return plot 

def indices_to_names_to_plot_tics(selection_indices:List[int], selection_dataset:dataset) -> hvplot.hvPlot:
  """helper function to bundle the selection's names as argument for plot_tics

  Parameters
  ----------
  selection_indices : List[int]
      selection of rows in the df_widget
  selection_dataset : dataset
      object of the dataset currently being displayed
  
  Returns
  -------
  hvplot.hvPlot
      hvPlot from calling plot_tics
  """
  names = selection_dataset.main_index_to_names(selection_indices)
  return plot_tics(names, selection_dataset.tics) 

def plot_runsticker(selection:List[int], data:pd.DataFrame) -> hv.Layout:
  """plots basic info about the selection of runs in a compact form

  Parameters
  ----------
  selection : List[int]
      list of indices(, selected from df_widget)
  data : pd.DataFrame
      main data frame containing the single value metrics (columns=[
        'Date', '# MS1', '# MS2', '# ID MS2', '# Peptidoforms', '# Proteoforms',
        '# Proteins', '# Features', '# ID Features', u'# Signal fluct. ↑',
        u'# Signal fluct. ↓', 'RT range left', 'RT range right', 
        'MZ range left', 'MZ range right', 'Name'])

  Returns
  -------
  hv.Layout
      a HoloViews layout of four plots:
      * idax = barplot('# MS1','# MS2', '# ID MS2')
      * qaax = barplot('# Features','# ID Features')
      * mzax = lineplot(m/z Range Setting)
      * flax = barplot(u'# Signal fluct. ↓', u'# Signal fluct. ↑')
  """
  qc2log.debug("plot_runsticker")
  if (data.shape[0] == 0) or\
     (len(selection)==0):
    blank = pd.DataFrame(columns=['Date','Counts', 'Ranges']).set_index('Date')\
      .hvplot.line(frame_height=200, frame_width=200).opts(xaxis='bare', yaxis='bare') * hv.Text(.5,.5, 'Select Runs \n from Table').opts(color='grey')
    return (blank + blank + blank + blank).cols(4).opts(shared_axes=False, toolbar='right', title='Counts & Ranges')

  qc2log.debug("plot_runsticker in ernest")
  # with pn.param.set_values(col[3][1], loading=True):
  while True:
    selection_df = data.iloc[selection].copy()
    selection_df.Date = selection_df.Date.dt.date
    idax = selection_df[['# MS1','# MS2', '# ID MS2', 'Date']].set_index('Date').hvplot.barh(
      frame_width=200, frame_height=200, xlabel='', ylabel="Count", title="Spectrum Counts")
    qaax = selection_df[['# Features','# ID Features', 'Date']].set_index('Date').hvplot.barh(
      frame_width=200, frame_height=200, xlabel='', ylabel="Count", title="Feature Counts")
    mzax = pd.melt(selection_df[['MZ range left', 'MZ range right', 'Date .. Name']].reset_index(), 
                  id_vars=['Date .. Name','index'], value_vars=['MZ range left','MZ range right']).hvplot.line(
                              x='value', y='index', by='Date .. Name', 
                              xlabel="m/z Range Setting", line_width=12,
                              frame_width=200, frame_height=200).opts(yaxis='bare', title="Acquisition Range")
    flax = selection_df[[u'# Signal fluct. ↓', u'# Signal fluct. ↑', 'Date']]\
            .set_index('Date').astype('int').hvplot.barh(
                frame_width=200, frame_height=200, xlabel='', ylabel="Count", title="Signal fluctuations")
    return (idax + qaax + mzax + flax).cols(4).opts(shared_axes=False, toolbar='right')

def plot_metrics_daterange(start:dt.datetime, end:dt.datetime, selection:List[str], data:pd.DataFrame, data_meta:pd.DataFrame) -> hvplot.hvPlot:
  """plots the selection of single value metrics over a selected daterange

  Parameters
  ----------
  start : dt.datetime
      start of the selected daterange
  end : dt.datetime
      end of the selected daterange 
  selection : List[str]
      list of column names representing the single value metrics 
  data : pd.DataFrame
      source data for the plot (columns=[
        'Date', '# MS1', '# MS2', '# ID MS2', '# Peptidoforms', '# Proteoforms',
        '# Proteins', '# Features', '# ID Features', u'# Signal fluct. ↑',
        u'# Signal fluct. ↓', 'RT range left', 'RT range right', 
        'MZ range left', 'MZ range right', 'Name'])
  data_meta: pd.DataFrame
      source data for the plot action annotations (columns={
                      "user_date": "Date", "user_email": "Contact", 
                      "additional_information": "Additional_Information", 
                      "problems": "Problems", "actions": "Actions"})

  Returns
  -------
  hvplot.hvPlot
      a lineplot overlay of all selected columns and selected dates
  """
  qc2log.debug("plot_metrics_daterange")
  if (data.shape[0] == 0) or\
     (len(selection)==0):
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(xlabel='Date', title='Single Metrics Timeseriesplot')

  qc2log.debug("plot_metrics_daterange in ernest")
  truncated = data.loc[:, (data.columns[
    (data.columns.str.startswith(('#','Date'))) & (data.columns.isin(selection+['Date']))
    ])]
  truncated.sort_values(by='Date', inplace = True)
  truncated = truncated[(truncated.Date.dt.date >= start.date()) & (truncated.Date.dt.date <= end.date())]
  if truncated.shape[0] == 0:
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(xlabel='Date', title='Single Metrics Timeseriesplot')

  truncated.Date = truncated.Date.dt.date
  line_plot = truncated.hvplot.line(x='Date', title='Single Metrics Timeseriesplot').opts(shared_axes=False)

  
  truncated_meta = data_meta[(data_meta.Date.dt.date >= start.date()) & (data_meta.Date.dt.date <= end.date())]
  annotations = [hv.VLine(a) for a in truncated_meta.Date.dt.date.to_list()]

  return hv.Overlay([line_plot] + annotations).opts(
    hv.opts.VLine(line_width=0.2, line_dash='dashed', color='#0096FF')
    ).opts(shared_axes=False, frame_height=400, frame_width=800)

def plot_ccharts(start:dt.datetime, end:dt.datetime, data_peps:pd.DataFrame, 
                 peps_mean: Dict[str,float], peps_std: Dict[str,float]) -> hv.Layout:
  """Xbar control charts for the mean deviation per run for ['dRT', 'dPPM', '# Chargestates', '# Identified QC2 Peptides']

  Parameters
  ----------
  start : dt.datetime
      start of the selected daterange
  end : dt.datetime
      end of the selected daterange 
  data_main : pd.DataFrame
      source data for the plot (columns=[
        'Date', '# MS1', '# MS2', '# ID MS2', '# Peptidoforms', '# Proteoforms',
        '# Proteins', '# Features', '# ID Features', u'# Signal fluct. ↑',
        u'# Signal fluct. ↓', 'RT range left', 'RT range right', 'MZ range left', 'MZ range right', 'Name'])
  data_peps : pd.DataFrame
      source data for the plot (columns=['dRT', 'dPPM', '# Chargestates', 
        '# Identified QC2 Peptides'])
  peps_mean: Dict[str,float]
      per dataset means of ['dRT', 'dPPM', '# Chargestates', '# Identified QC2 Peptides'] packed in a dict
  peps_std: Dict[str,float]
      per dataset standard deviations of ['dRT', 'dPPM', '# Chargestates', '# Identified QC2 Peptides'] packed in a dict

  Returns
  -------
  hv.Layout
      a HoloViews layout of 4 control charts
  """
  qc2log.debug("plot_ccharts")
  if (data_peps.shape[0] == 0) or\
     (peps_mean is None) or\
     (peps_std is None):
    layout = hv.Layout([pd.DataFrame(columns=["Date","mean per day"]).hvplot.line()]*4).cols(2).opts(shared_axes=False)
    return layout
    
  qc2log.debug("plot_metrics_daterange in ernest")
  pep_df = data_peps
  filtered_pep_df = pep_df[(pep_df['Date'].dt.date >=  start.date()) & (pep_df['Date'].dt.date <=  end.date())]
  filtered_pep_df_means = pd.DataFrame({
      "dRT": filtered_pep_df.groupby(['Date','Name'])['dRT'].mean(),
      "dPPM": filtered_pep_df.groupby(['Date','Name'])['dPPM'].mean(),
      "# Chargestates": filtered_pep_df.groupby(['Date','Name'])['Chargestate'].nunique(),
      "# Identified QC2 Peptides": filtered_pep_df.groupby(['Date','Name'])['Peptide'].nunique(),
  })

  figs = list()
  for metric_name in filtered_pep_df_means.keys():
    hline = hv.HLine(peps_mean[metric_name])
    hline.opts(
        color='red',
        line_dash='dashed',
        line_width=2.0,
    )

    hspans_t0 = hv.HSpan(peps_mean[metric_name], peps_mean[metric_name]+1*peps_std[metric_name]).opts(color='#cadeab')
    hspans_b0 = hv.HSpan(peps_mean[metric_name], peps_mean[metric_name]-1*peps_std[metric_name]).opts(color='#cadeab')
    hspans_t1 = hv.HSpan(peps_mean[metric_name]+1*peps_std[metric_name], peps_mean[metric_name]+2*peps_std[metric_name]).opts(color='#ffd699')
    hspans_b1 = hv.HSpan(peps_mean[metric_name]-1*peps_std[metric_name], peps_mean[metric_name]-2*peps_std[metric_name]).opts(color='#ffd699')
    hspans_t2 = hv.HSpan(peps_mean[metric_name]+2*peps_std[metric_name], peps_mean[metric_name]+3*peps_std[metric_name]).opts(color='#933126')
    hspans_b2 = hv.HSpan(peps_mean[metric_name]-2*peps_std[metric_name], peps_mean[metric_name]-3*peps_std[metric_name]).opts(color='#933126')

    fig = hspans_t0 * hspans_t1 * hspans_t2 * hspans_b0 * hspans_b1 * hspans_b2 * hline * \
      filtered_pep_df_means.hvplot.line(y=metric_name, title= metric_name + " per day mean cchart")
    fig.opts(ylim=(peps_mean[metric_name]-3*peps_std[metric_name], peps_mean[metric_name]+3*peps_std[metric_name]),
             default_tools=[], active_tools=[], tools=['wheel_zoom', 'save', 'reset', 'hover'])
    figs.append(fig)
  layout = hv.Layout(figs).cols(2).opts(shared_axes=False)
  return layout

class dataset_panels:
  """generates the panels and widgets for a given dataset and holds them"""
  DF_MAX_SELECT = 6
  def __init__(self, dataset: dataset):
    self.ds_descriptor = pn.pane.Markdown(name="Dataset descriptor", 
                                          object=DESCR_TEMPL.format(s=dataset.start_str(),
                                            e=dataset.end_str(),
                                            p=dataset.label,
                                            n=dataset.size_str(),
                                            i=dataset.instrument))

    self.df_widget = pn.widgets.Tabulator(value=dataset.main,
                                    selectable=dataset_panels.DF_MAX_SELECT,
                                    show_index=False,
                                    pagination='local', page_size=10,
                                    disabled=True,
                                    hidden_columns=["Date", "Name",
                                                    'RT range left', 'RT range right',
                                                    'MZ range left','MZ range right'])

    self.date_range_slider = pn.widgets.DateRangeSlider(
      name='Date Range Slider',
      start=dt.datetime.combine(dataset.main.Date.min(),dt.datetime.min.time()),
      end=dt.datetime.combine(dataset.main.Date.max(),dt.datetime.min.time()),
      # value=( dt.datetime.now()-dt.timedelta(days=1), dt.datetime.now()+dt.timedelta(days=1) ),
      tooltips=True)

    self.checkbox_group = pn.widgets.CheckBoxGroup(
      options=dataset.main.columns.drop(['Date', 'Name','Date .. Name', 
                                         'RT range left','RT range right', 
                                         'MZ range left','MZ range right']).to_list(),
      value=['# MS1', '# MS2'], inline=False, name='Checkbox Group', )
    
    self.calendar_hist_pane = pn.pane.HoloViews(name="Calendar histogram pane",
                                        object=plot_calendar_hist(dataset.main))

    self.tics_pane = pn.bind(indices_to_names_to_plot_tics, 
                      selection_indices=self.df_widget.param.selection,  
                      selection_dataset=dataset)

    self.runsticker_pane = pn.bind(plot_runsticker,
                            selection=self.df_widget.param.selection, 
                            data=dataset.main)

    self.metrics_pane = pn.bind(plot_metrics_daterange, 
                          start=self.date_range_slider.param.value_start, 
                          end=self.date_range_slider.param.value_end,
                          selection=self.checkbox_group.param.value, 
                          data=dataset.main,data_meta=dataset.meta)
    
    self.cchart_pane = pn.bind(plot_ccharts, 
                        start=self.date_range_slider.param.value_start, 
                        end=self.date_range_slider.param.value_end,
                        data_peps=dataset.peps,
                        peps_mean=dataset.peps_mean, 
                        peps_std=dataset.peps_std)

# """main loop code starts here"""

# mzqc_basepath = "/content/drive/Shareddrives/mzqclib-manuscript/test data/CRG"
mzqc_basepath = "mzqcs"
# mzqc_basepath = "smalldevset"
dataset_mid_path = {"Lumos_2017": "PXD019888",
"Velos_2018": "PXD019889",
"Lumos_2018": "PXD019891",
"Velos_2017": "PXD019892"}
metadata_paths = {k: os.path.join(os.path.join(mzqc_basepath, v), "metadata_"+k+".json") for k,v in dataset_mid_path.items()}
mzqc_paths = {k: [os.path.join(os.path.join(mzqc_basepath, v), x) for x in os.listdir(os.path.join(mzqc_basepath, v)) if x.endswith(".mzqc")] for k,v in dataset_mid_path.items()}

ds_select = pn.widgets.Select(name='Select Dataset', options=list(mzqc_paths.keys()))
app_col = None
def update_ds(ds_key: str, mzqc_paths: Dict[str,str]):
  global app_col
  global ds_switch
  if app_col:
    with pn.param.set_values(app_col[1], loading=True),\
        pn.param.set_values(app_col[2], loading=True),\
        pn.param.set_values(app_col[3], loading=True),\
        pn.param.set_values(app_col[4], loading=True),\
        pn.param.set_values(app_col[5], loading=True),\
        pn.param.set_values(app_col[6], loading=True):
      current_dataset = load_ds(ds_key=ds_key, mzqc_paths=mzqc_paths, metadata_paths=metadata_paths)
      current_panels = dataset_panels(current_dataset)
      row1 = pn.Row(current_panels.ds_descriptor,
                    pn.Spacer(width=100, height=300), 
                    current_panels.calendar_hist_pane, align='center')
      row2 = pn.Row(current_panels.df_widget)
      row3 = pn.Row(current_panels.tics_pane)
      row4 = pn.Row(pn.Spacer(width=20, height=300), 
                    current_panels.runsticker_pane)
      internal_col = pn.Column(current_panels.date_range_slider, 
                              current_panels.checkbox_group)
      row5 = pn.Row(pn.Spacer(width=150, height=300), 
                    internal_col, 
                    current_panels.metrics_pane)
      row6 = pn.Row(pn.Tabs(('Xbar chart',current_panels.cchart_pane), 
                            ('Rbar chart',pn.Row(current_panels.cchart_pane)),  
                            tabs_location='left'))   # TODO ,('Rbar chart', !cchart_pane)
      app_col[1:]=[row1,row2,row3,row4,row5,row6]

ds_switch = pn.bind(update_ds, 
                    ds_key=ds_select.param.value, 
                    mzqc_paths=mzqc_paths)

app_col = pn.Column(pn.Row(pn.pane.Markdown(name="Dataset descriptor", object="# QC2 View"),
                           ds_select,ds_switch), 
                    pn.Row(pn.Spacer(styles=dict(background='WhiteSmoke'), sizing_mode='stretch_both'), height=200, width=1000), 
                    pn.Row(pn.Spacer(styles=dict(background='WhiteSmoke'), sizing_mode='stretch_both'), height=200, width=1000), 
                    pn.Row(pn.Spacer(styles=dict(background='WhiteSmoke'), sizing_mode='stretch_both'), height=200, width=1000), 
                    pn.Row(pn.Spacer(styles=dict(background='WhiteSmoke'), sizing_mode='stretch_both'), height=200, width=1000), 
                    pn.Row(pn.Spacer(styles=dict(background='WhiteSmoke'), sizing_mode='stretch_both'), height=200, width=1000), 
                    pn.Row(pn.Spacer(styles=dict(background='WhiteSmoke'), sizing_mode='stretch_both'), height=200, width=1000))
app_col.servable(title="QC2 Dashboard")

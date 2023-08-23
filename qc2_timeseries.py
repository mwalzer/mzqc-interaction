# -*- coding: utf-8 -*-
"""CRG-QC2-timeseries.ipynb

## Load Dependencies and Data
"""
import os
from typing import Dict, List, Tuple

import pandas as pd
from mzqc import MZQCFile as qc
import datetime as dt

import panel as pn
import matplotlib.pyplot as plt
import hvplot.pandas
import holoviews as hv
pn.extension('tabulator')
hv.extension('bokeh')

from dataclasses import dataclass, field

@dataclass
class dataset:
    """Class for keeping track the different metric dataframes for a particular dataset."""
    main: pd.DataFrame = pd.DataFrame(columns=['Date','# MS1','# MS2', '# ID MS2','# Peptidoforms',
                                '# Proteoforms','# Proteins','# Features',
                                '# ID Features',u'# Signal fluct. ↑',
                                u'# Signal fluct. ↓','RT range','MZ range','mzrange',
                                'Name'])  # closely track the load function
    tics: pd.DataFrame = pd.DataFrame(columns=["RT", "Intensity", "# Peak", "NativeID", "Date", "Name"])  # closely track the load function 
    peps: pd.DataFrame = pd.DataFrame(columns=["RT", "Peptide", "dPPM", "Chargestate", "Date", "Name"])  # closely track the load function
    peps_std: Dict[str,float] = field(default_factory=dict)
    peps_mean: Dict[str,float] = field(default_factory=dict)
    label: str = ""
    instrument: str = ""

    def start_str(self):
      return str(self.main.Date.min().date())
    def end_str(self):
      return str(self.main.Date.max().date())
    def size_str(self):
      return str(len(self.main))


# mzqc_basepath = "/content/drive/Shareddrives/mzqclib-manuscript/test data/CRG"
# mzqc_basepath = "mzqcs"
mzqc_basepath = "smalldevset"
dataset_mid_path = {"Lumos_2017": "PXD019888",
"Velos_2018": "PXD019889",
"Lumos_2018": "PXD019891",
"Velos_2017": "PXD019892"}
metadata_paths = {k: os.path.join(os.path.join(mzqc_basepath, v), "metadata_"+k+".json") for k,v in dataset_mid_path.items()}
mzqc_paths = {k: [os.path.join(os.path.join(mzqc_basepath, v), x) for x in os.listdir(os.path.join(mzqc_basepath, v)) if x.endswith(".mzqc")] for k,v in dataset_mid_path.items()}

df_widget = None
date_range_slider = None
checkbox_group = None 
ds_descriptor = None
current_dataset = dataset()

def extract_date(rowrunmzqc: qc.MzQcFile):
  rawdate = next(iter([pd.to_datetime(cd.value) for cd in rowrunmzqc.runQualities[0].metadata.inputFiles[0].fileProperties if cd.accession=="MS:1000747"]))
  return pd.Timestamp(rawdate)

def load_main_from_mzqc(rowrunname: str, rowrunmzqc: qc.MzQcFile):
  return {'Date': extract_date(rowrunmzqc),
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
          'RT range': str(next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000070"]))),
          'MZ range': str(next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000069"]))),
          'mzrange': next(iter([cd.value for cd in rowrunmzqc.runQualities[0].qualityMetrics if cd.accession=="MS:4000069"])),
          'Name':rowrunname,}

def load_tic_from_mzqc(rowrunname: str, rowrunmzqc: qc.MzQcFile):
  date = extract_date(rowrunmzqc)
  tics = pd.DataFrame(next(iter([m.value for m in rowrunmzqc.runQualities[0].qualityMetrics if m.accession=="MS:4000104"])))\
              .rename(columns={"MS:1000894": "RT", "MS:1000285": "Intensity", "MS:1003059": "# Peak", "MS:1000767": "NativeID"})\
              .assign(**{"Date":date, "Name":rowrunname})
  return tics
  
def load_peps_from_mzqc(rowrunname: str, rowrunmzqc: qc.MzQcFile):
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

def load_ds(ds_key: str, mzqc_paths: Dict[str,str], dataset: dataset, ):
  if ds_key not in mzqc_paths.keys():
    print("ds_key not in mzqc_paths.")
    return
  print(">>>updating with {}".format(ds_key))

  dfs_main = list()
  dfs_tics = list()
  dfs_peps = list()
  for mzqc_path in mzqc_paths.get(ds_key, None):
    with open(mzqc_path, "r") as file:
        mzqcobj = qc.JsonSerialisable.FromJson(file)
        name = os.path.basename(mzqc_path)
        dfs_main.append(load_main_from_mzqc(name, mzqcobj))
        dfs_tics.append(load_tic_from_mzqc(name, mzqcobj))
        dfs_peps.append(load_peps_from_mzqc(name,mzqcobj))
  dataset.main = pd.DataFrame(dfs_main)
  dataset.peps = pd.concat(dfs_peps, axis=0)
  dataset.tics = pd.concat(dfs_tics, axis=0)
  dataset.peps = pd.merge(dataset.peps, pd.DataFrame(
    {"mean_rt_pp": dataset.peps.groupby(['Peptide'])['RT'].mean()}), how="inner", on="Peptide")
  dataset.peps["dRT"] = dataset.peps['mean_rt_pp'] - dataset.peps['RT'].drop(columns=["mean_rt_pp"])
  # ds_main = ds_main.astype({'Date':'datetime64[ns]'})
  dataset.main.sort_values(by='Date', inplace = True) 
  dataset.peps.sort_values(by='Date', inplace = True) 
  dataset.tics.sort_values(by='Date', inplace = True) 

  dataset.peps_std = {
      "dRT": dataset.peps['dRT'].std(),
      "dPPM": dataset.peps['dPPM'].std(),
      "# Chargestates": dataset.peps.groupby(['Date','Name'])['Chargestate'].nunique().std(),
      "# Identified QC2 Peptides": dataset.peps.groupby(['Date','Name'])['Peptide'].nunique().std(),
  }
  dataset.peps_mean = {
      "dRT": dataset.peps['dRT'].mean(),
      "dPPM": dataset.peps['dPPM'].mean(),
      "# Chargestates": dataset.peps.groupby(['Date','Name'])['Chargestate'].nunique().mean(),
      "# Identified QC2 Peptides": dataset.peps.groupby(['Date','Name'])['Peptide'].nunique().mean(),
  }

  dataset.label = dataset_mid_path[ds_key]
  dataset.instrument = "Thermo {}".format(next(iter(ds_key.split('_'))))

  update_from_ds_load(current_dataset)
  return

def plot_tics(selection:List[int], update_key:bool, dataset:dataset):
  print(">>>plot_tics")
  if (dataset is None) or\
     (len(selection)==0) or\
     (dataset.main.shape[0]==0):
    return pd.DataFrame(columns=["RT", "Intensity"]).hvplot.line(title="Total Ion Chromatogram",
                                                            line_width=0.5, 
                                                            xlabel="Retentiontime [s]", 
                                                            ylabel="Relative Intensity of Ion Current",
                                                            frame_height=500,
                                                            frame_width=800,)
  
  print(">>>plot_tics in ernest")
  # with pn.param.set_values(col[3][0], loading=True):
  while True:
    ds_tic = dataset.tics
    selected_tic = dataset.tics[ds_tic.Name.isin(dataset.main.iloc[selection].Name)]
    selected_tic["str_date"] = selected_tic.Date.dt.strftime('%Y-%m-%d')
    selected_tic["Date .. Name"] = selected_tic.str_date + ".." + selected_tic.Name.str[-10:]
    selected_tic["Date .. Name"] = selected_tic["Date .. Name"].astype(str)
    selected_tic = selected_tic[["RT","Intensity","Date .. Name"]].sort_values(["RT"], axis = 0, ascending = True)
    plot = selected_tic.hvplot.line(x="RT",y="Intensity", groupby=["Date .. Name"],
                                    title="Total Ion Chromatogram",
                                    xlabel="Retentiontime [s]", ylabel="Relative Intensity of Ion Current",
                                    frame_height=500, frame_width=800,).overlay()
    return plot 

def plot_runsticker(selection:List[int], update_key:bool, dataset:dataset):
  print(">>>plot_runsticker")
  if (dataset is None) or\
      (len(selection)==0) or\
      (dataset.main.shape[0] == 0):
    selection_df = pd.DataFrame(columns=['Date','# MS1','# MS2','# ID MS2','# Features','# ID Features', '# Signal fluct. ↓', '# Signal fluct. ↑',])
    idax = selection_df[['# MS1','# MS2', '# ID MS2', 'Date']].set_index('Date').hvplot.barh(frame_height=200, frame_width=200)
    qaax = selection_df[['# Features','# ID Features', 'Date']].set_index('Date').hvplot.barh(frame_height=200, frame_width=200)
    mzax = selection_df.hvplot.line(x="Date", xlabel="m/z Range Setting", ylim=(100,1600),
                             invert=True, frame_height=200, frame_width=200).opts(yaxis='bare')
    flax = selection_df[['# Signal fluct. ↓', '# Signal fluct. ↑', 'Date']]\
            .set_index('Date').astype('int').hvplot.bar(rot=45, frame_height=200, frame_width=200)
    return (idax + qaax + mzax + flax).cols(2).opts(shared_axes=False)

  print(">>>plot_runsticker in ernest")
  # with pn.param.set_values(col[3][1], loading=True):
  while True:
    selection_df = dataset.main.iloc[selection]
    idax = selection_df[['# MS1','# MS2', '# ID MS2', 'Date']].set_index('Date').hvplot.barh(frame_width=200, frame_height=200)
    qaax = selection_df[['# Features','# ID Features', 'Date']].set_index('Date').hvplot.barh(frame_width=200, frame_height=200)
    tmpdf = pd.pivot(selection_df[['mzrange', 'Date']]\
          .reset_index().rename(columns={'index':'xpos'})\
          .explode('mzrange', ignore_index=True).reset_index(), index=["index","xpos"]\
          , columns="Date", values="mzrange").reset_index(level=("xpos",))
    mzax = tmpdf.hvplot.line(x="xpos",  
                             xlim=(tmpdf.xpos.min()-0.5,tmpdf.xpos.max()+0.5),
                             xlabel="m/z Range Setting",
                             ylim=(100,1600),
                             invert=True, 
                             frame_width=200,
                             frame_height=200).opts(yaxis='bare')
    flax = selection_df[['# Signal fluct. ↓', '# Signal fluct. ↑', 'Date']]\
            .set_index('Date').astype('int').hvplot.bar(rot=45, frame_width=200, frame_height=200)
    return (idax + qaax + mzax + flax).cols(2).opts(shared_axes=False,)

def plot_metrics_daterange(start:dt.datetime, end:dt.datetime, selection:List[str], update_key:bool, dataset:dataset):
  print(">>>plot_metrics_daterange")
  if (dataset is None) or\
     (len(selection)==0) or\
     (dataset.main.shape[0] == 0):
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(xlabel='Date', title='Single Metrics Timeseriesplot')

  print(">>>plot_metrics_daterange in ernest")
  truncated = dataset.main.loc[:, (dataset.main.columns[
    (dataset.main.columns.str.startswith(('#','Date'))) & (dataset.main.columns.isin(selection+['Date']))
    ])]
  truncated.sort_values(by='Date', inplace = True)
  print("df", type(truncated.Date[0]))  
  print("slider", type(start),start,end)
  truncated = truncated[(truncated.Date.dt.date >= start.date()) & (truncated.Date.dt.date <= end.date())]
  print(truncated)
  if truncated.shape[0] == 0:
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(xlabel='Date', title='Single Metrics Timeseriesplot')

  truncated.Date = truncated.Date.dt.date
  line_plot = truncated.hvplot.line(x='Date', title='Single Metrics Timeseriesplot')
  print("boingboing")
  return line_plot

def plot_ccharts(start, end, update_key:bool, dataset:dataset):
  if (dataset is None) or\
     (dataset.main.shape[0] == 0) or\
     (any([getattr(dataset,x,None).shape[0] ==0 for x in ['peps','peps_mean','peps_std']])):
    layout = hv.Layout([pd.DataFrame(columns=["Date","mean per day"]).hvplot.line()]*4).cols(2).opts(shared_axes=False)
    return layout
  
  pep_df = dataset.peps

  filtered_pep_df = pep_df[(pep_df['Date'].dt.date >=  start.date()) & (pep_df['Date'].dt.date <=  end.date())]
  filtered_pep_df_means = pd.DataFrame({
      "dRT": filtered_pep_df.groupby(['Date','Name'])['dRT'].mean(),
      "dPPM": filtered_pep_df.groupby(['Date','Name'])['dPPM'].mean(),
      "# Chargestates": filtered_pep_df.groupby(['Date','Name'])['Chargestate'].nunique(),
      "# Identified QC2 Peptides": filtered_pep_df.groupby(['Date','Name'])['Peptide'].nunique(),
  })

  figs = list()
  pep_df_mean = dataset.peps_mean
  pep_df_std = dataset.peps_std
  for metric_name in filtered_pep_df_means.keys():
    hline = hv.HLine(pep_df_mean[metric_name])
    hline.opts(
        color='red',
        line_dash='dashed',
        line_width=2.0,
    )

    hspans_t0 = hv.HSpan(pep_df_mean[metric_name], pep_df_mean[metric_name]+1*pep_df_std[metric_name]).opts(color='#cadeab')
    hspans_b0 = hv.HSpan(pep_df_mean[metric_name], pep_df_mean[metric_name]-1*pep_df_std[metric_name]).opts(color='#cadeab')
    hspans_t1 = hv.HSpan(pep_df_mean[metric_name]+1*pep_df_std[metric_name], pep_df_mean[metric_name]+2*pep_df_std[metric_name]).opts(color='#ffd699')
    hspans_b1 = hv.HSpan(pep_df_mean[metric_name]-1*pep_df_std[metric_name], pep_df_mean[metric_name]-2*pep_df_std[metric_name]).opts(color='#ffd699')
    hspans_t2 = hv.HSpan(pep_df_mean[metric_name]+2*pep_df_std[metric_name], pep_df_mean[metric_name]+3*pep_df_std[metric_name]).opts(color='#933126')
    hspans_b2 = hv.HSpan(pep_df_mean[metric_name]-2*pep_df_std[metric_name], pep_df_mean[metric_name]-3*pep_df_std[metric_name]).opts(color='#933126')

    fig = hspans_t0 * hspans_t1 * hspans_t2 * hspans_b0 * hspans_b1 * hspans_b2 * hline * \
      filtered_pep_df_means.hvplot.line(y=metric_name, title= metric_name + " per day mean cchart")
    fig.opts(ylim=(pep_df_mean[metric_name]-3*pep_df_std[metric_name], pep_df_mean[metric_name]+3*pep_df_std[metric_name]),
             default_tools=[], active_tools=[], tools=['wheel_zoom', 'save', 'reset', 'hover'],toolbar='above')
    figs.append(fig)

  layout = hv.Layout(figs).cols(2).opts(shared_axes=False)
  return layout

def plot_calendar_hist(update_key:bool, dataset:dataset):
  if (dataset is None) or\
     (dataset.main.shape[0] == 0):
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(
          title="QC2 Run Yearly Distribution", 
          ylabel="# Runs",
          xlabel="Month of {}".format("N/A"))

  calendar_hist = dataset.main.assign(month=dataset.main.Date.dt.month)\
        .groupby(by=['month']).count().hvplot.bar(y="Date")\
        .opts(title="QC2 Run Yearly Distribution", ylabel="# Runs",
              xlabel="Month of {}".format(next(iter(dataset.main.Date.dt.year))))
  return calendar_hist

def update_from_ds_load(dataset:dataset):
  global df_widget
  global ds_descriptor
  global date_range_slider
  global checkbox_group
  global update_ds_triggered
  print(">>>updating depended")

  df_widget.value = dataset.main
  df_widget.selection = [0]
  print(">>>updated df_widget")

  ds_descriptor.object = DESCR_TEMPL.format(s=dataset.start_str(),
                                          e=dataset.end_str(),
                                          p=dataset.label,
                                          n=dataset.size_str(),
                                          i=dataset.instrument,)
  print(">>>updated dataset description")

  # update date_range_slider 
  print(">>> pre update_ds date_range_slider types", 
        type(date_range_slider.start), 
        type(date_range_slider.end),
        date_range_slider.start, 
        date_range_slider.end)
  date_range_slider.start = dt.datetime.combine(dataset.main.Date.min(),dt.datetime.min.time())
  date_range_slider.end = dt.datetime.combine(dataset.main.Date.max(),dt.datetime.min.time())
  print(">>> after update_ds date_range_slider types", 
        type(date_range_slider.start), 
        type(date_range_slider.end),
        date_range_slider.start, 
        date_range_slider.end)
  date_range_slider.value = (date_range_slider.start,date_range_slider.end)  
  # the last line causes trouble because it triggers the bind to plot_metrics_daterange
  # maybe if the value set is done after 1.switching main 2.updating to default plot 3.set date_range_slider values
  # or maybe use jslink https://panel.holoviz.org/how_to/links/link_plots.html
  print(">>>updated date_range_slider")

  # update checkbox_group
  setattr(checkbox_group, 'options', 
          dataset.main.columns.drop(['Date', 'RT range', 'MZ range', 'mzrange', 'Name']).to_list())  
          # do not deactivate unless all mzQC produce the same column layout
  setattr(checkbox_group, 'value', ['# MS1', '# MS2'])
  print(">>>updated checkbox_group")

  # update_ds_triggered.value = not update_ds_triggered.value  # this to trigger calendar_hist_pane will break panel, also metrics_pane does not need that to update?

  print(">>>updated widgets after dataset load!")
  return

df_widget_max_select = 6
df_widget = pn.widgets.Tabulator(value=current_dataset.main,
                                  selectable='checkbox',
                                  show_index=False,
                                  pagination='local', page_size=10,
                                  disabled=True,
                                  hidden_columns=["Name",'mzrange', 'RT Range', 'MZ Range'])
df_widget.selectable=df_widget_max_select

date_range_slider = pn.widgets.DateRangeSlider(
    name='Date Range Slider',
    start=dt.datetime.now()-dt.timedelta(days=1), end=dt.datetime.now()+dt.timedelta(days=1),
    value=( dt.datetime.now()-dt.timedelta(days=1), dt.datetime.now()+dt.timedelta(days=1) ),
    tooltips=True)

checkbox_group = pn.widgets.CheckBoxGroup(
    name='Checkbox Group', value=[], 
    options=[],
    inline=False
)

DESCR_TEMPL = """
# The Data Set

This dataset ({p}) is comprised of mzQC files from {n} runs,

which were measured over the timecourse between {s} and {e}, 

from QC2 samples on an instrument of type {i}.

You can in inspect the details in the following interactive charts.
"""
ds_descriptor = pn.pane.Markdown(DESCR_TEMPL.format(s="1900-01-01",
                                                    e="1900-01-01",
                                                    p="PXD------",
                                                    n=str(len([])),
                                                    i="\<Make\> \<Model\>"))

ds_select = pn.widgets.Select(name='Select Dataset', options=list(mzqc_paths.keys()))
update_ds_triggered = pn.widgets.Toggle(name="update_ds_triggered")

plot_tics
tics_pane = pn.bind(plot_tics, 
                    selection=df_widget.param.selection, 
                    update_key=update_ds_triggered.param.value, 
                    dataset=current_dataset)
runsticker_pane = pn.bind(plot_runsticker,
                          selection=df_widget.param.selection, 
                       update_key=update_ds_triggered.param.value, 
                       dataset=current_dataset)

metrics_pane = pn.bind(plot_metrics_daterange, 
                       start=date_range_slider.param.value_start, 
                       end=date_range_slider.param.value_end,
                       selection=checkbox_group.param.value, 
                       update_key=update_ds_triggered.param.value, 
                       dataset=current_dataset)
cchart_pane = pn.bind(plot_ccharts, 
                      start=date_range_slider.param.value_start, 
                      end=date_range_slider.param.value_end,
                      update_key=update_ds_triggered.param.value, 
                      dataset=current_dataset)

calendar_hist_pane = pn.bind(plot_calendar_hist, 
                             update_key=update_ds_triggered.param.value, 
                             dataset=current_dataset,)

ds_switch = pn.bind(load_ds, 
                    ds_key=ds_select.param.value, 
                    mzqc_paths=mzqc_paths,
                    dataset=current_dataset, 
                    )

row0 = pn.Row(ds_select)
row1 = pn.Row(ds_descriptor, calendar_hist_pane)
row2 = pn.Row(df_widget)
row3 = pn.Row(tics_pane, runsticker_pane)
internal_col = pn.Column(date_range_slider, checkbox_group)
row4 = pn.Row(internal_col, metrics_pane)
# row5 = pn.Row(tabs = pn.Tabs(cchart_pane, tabs_location='left'))   #,('Rbar control chart', cchart_pane)

col = pn.Column(ds_switch) # bind must be included to be active
col.append(row0)  
col.append(row1) 
col.append(row2)
col.append(row3)
col.append(row4)
# col.append(row5)
# col.append(cchart_pane)
col.servable()



# ! python /home/vscode/.local/bin/panel serve crg_qc2_timeseries.py
# ! panel serve qc2_timeseries.py
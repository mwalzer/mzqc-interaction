# -*- coding: utf-8 -*-
"""CRG-QC2-timeseries.ipynb

## Load Dependencies and Data
"""
import os
import re
from typing import Dict, List

import pandas as pd
from mzqc import MZQCFile as qc
import datetime as dt

import panel as pn
import matplotlib.pyplot as plt
import hvplot.pandas
import holoviews as hv
pn.extension('tabulator')
hv.extension('bokeh')

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

def load_ds(ds_key, mzqc_paths=mzqc_paths):
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

  ds_main = pd.DataFrame(dfs_main)
  ds_pep = pd.concat(dfs_peps, axis=0)
  ds_tic = pd.concat(dfs_tics, axis=0)
  # ds_main = ds_main.astype({'Date':'datetime64[ns]'})
  ds_main.sort_values(by='Date', inplace = True) 
  ds_pep.sort_values(by='Date', inplace = True) 
  ds_tic.sort_values(by='Date', inplace = True) 
  
  return ds_main, ds_tic, ds_pep

def plot_tics(selection=[], df_widget:pn.widgets.Tabulator=df_widget):
  if (df_widget is None) or (len(selection)==0) or (type(df_widget.value) != pd.DataFrame) or (df_widget.value.shape[0]):
    return pd.DataFrame(columns=["RT", "Intensity"]).hvplot.line(title="Total Ion Chromatogram",
                                                            line_width=0.5, 
                                                            xlabel="Retentiontime [s]", 
                                                            ylabel="Relative Intensity of Ion Current",
                                                            frame_height=500,
                                                            frame_width=800,)
  
  with pn.param.set_values(col[3][0], loading=True):
    ds_tic = df_widget.tics
    selected_tic = ds_tic[ds_tic.Name.isin(df_widget.value.iloc[selection].Name)]
    selected_tic["str_date"] = selected_tic.Date.dt.strftime('%Y-%m-%d')
    selected_tic["Date .. Name"] = selected_tic.str_date + ".." + selected_tic.Name.str[-10:]
    selected_tic["Date .. Name"] = selected_tic["Date .. Name"].astype(str)
    selected_tic = selected_tic[["RT","Intensity","Date .. Name"]].sort_values(["RT"], axis = 0, ascending = True)
    plot = selected_tic.hvplot.line(x="RT",y="Intensity", groupby=["Date .. Name"],
                                    title="Total Ion Chromatogram",
                                    xlabel="Retentiontime [s]", ylabel="Relative Intensity of Ion Current",
                                    frame_height=500, frame_width=800,).overlay()
    return plot 

def plot_runsticker(selection=[], df_widget:pn.widgets.Tabulator=df_widget):
  if (df_widget is None) or (len(selection)==0) or (type(df_widget.value) != pd.DataFrame) or (df_widget.value.shape[0]):
    selection_df = pd.DataFrame(columns=['Date','# MS1','# MS2','# ID MS2','# Features','# ID Features', '# Signal fluct. ↓', '# Signal fluct. ↑',])
    idax = selection_df[['# MS1','# MS2', '# ID MS2', 'Date']].set_index('Date').hvplot.barh(frame_height=200, frame_width=200)
    qaax = selection_df[['# Features','# ID Features', 'Date']].set_index('Date').hvplot.barh(frame_height=200, frame_width=200)
    mzax = selection_df.hvplot.line(x="Date", xlabel="m/z Range Setting", ylim=(100,1600),
                             invert=True, frame_height=200, frame_width=200).opts(yaxis='bare')
    flax = selection_df[['# Signal fluct. ↓', '# Signal fluct. ↑', 'Date']]\
            .set_index('Date').astype('int').hvplot.bar(rot=45, frame_height=200, frame_width=200)
    return (idax + qaax + mzax + flax).cols(2).opts(shared_axes=False)

  with pn.param.set_values(col[3][1], loading=True):
    selection_df = df_widget.value.iloc[selection]
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

def plot_metrics_daterange(start, end, selection, df_widget:pn.widgets.Tabulator=df_widget):
  if (df_widget is None) or (len(selection)==0) or (type(df_widget.value) != pd.DataFrame) or (df_widget.value.shape[0]):
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(xlabel='Date', title='Single Metrics Timeseriesplot')

  df = df_widget.value
  truncated = df.loc[:, (df.columns[(df.columns.str.startswith(('#','Date'))) & (df.columns.isin(selection+['Date']))])]
  truncated.sort_values(by='Date', inplace = True)
  print("df", type(truncated.Date[0]))  
  print("slider", type(start))
  truncated = truncated[(truncated.Date.dt.date >= start.date()) & (truncated.Date.dt.date <= end.date())]
  print(truncated)
  line_plot = truncated.hvplot.line(x='Date', title='Single Metrics Timeseriesplot')
  return line_plot

def plot_ccharts(start, end, df_widget:pn.widgets.Tabulator=df_widget):
  if (df_widget is None) or (type(df_widget.value) != pd.DataFrame) or (df_widget.value.shape[0]==0) or\
    (any([getattr(df_widget,x,None) is None for x in ['peps','peps_mean','peps_std']])):
    layout = hv.Layout([pd.DataFrame(columns=["Date","mean per day"]).hvplot.line()]*4).cols(2).opts(shared_axes=False)
    return layout
  
  pep_df = df_widget.peps

  filtered_pep_df = pep_df[(pep_df['Date'].dt.date >=  start.date()) & (pep_df['Date'].dt.date <=  end.date())]
  filtered_pep_df_means = pd.DataFrame({
      "dRT": filtered_pep_df.groupby(['Date','Name'])['dRT'].mean(),
      "dPPM": filtered_pep_df.groupby(['Date','Name'])['dPPM'].mean(),
      "# Chargestates": filtered_pep_df.groupby(['Date','Name'])['Chargestate'].nunique(),
      "# Identified QC2 Peptides": filtered_pep_df.groupby(['Date','Name'])['Peptide'].nunique(),
  })

  figs = list()
  pep_df_mean = df_widget.peps_mean
  pep_df_std = df_widget.peps_std
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

    fig = hspans_t0 * hspans_t1 * hspans_t2 * hspans_b0 * hspans_b1 * hspans_b2 * filtered_pep_df_means.hvplot.line(y=metric_name, title= metric_name + " per day mean cchart") * hline
    fig.opts(ylim=(pep_df_mean[metric_name]-3*pep_df_std[metric_name], pep_df_mean[metric_name]+3*pep_df_std[metric_name]),
             default_tools=[], active_tools=[], tools=['wheel_zoom', 'save', 'reset', 'hover'],toolbar='above')

    figs.append(fig)
  layout = hv.Layout(figs).cols(2).opts(shared_axes=False)
  return layout

def update_ds_pep_stats(df_widget:pn.widgets.Tabulator=df_widget):
  pep_df = df_widget.peps
  tdf = pd.merge(pep_df, pd.DataFrame({"mean_rt_pp": pep_df.groupby(['Peptide'])['RT'].mean()}), how="inner", on="Peptide")
  tdf["dRT"] = tdf['mean_rt_pp'] - tdf['RT'].drop(columns=["mean_rt_pp"])
  setattr(df_widget, 'peps', tdf)

  pep_df_std = {
      "dRT": tdf['dRT'].std(),
      "dPPM": tdf.dPPM.std(),
      "# Chargestates": tdf.groupby(['Date','Name'])['Chargestate'].nunique().std(),
      "# Identified QC2 Peptides": tdf.groupby(['Date','Name'])['Peptide'].nunique().std(),
  }
  setattr(df_widget, 'peps_std', pep_df_std)

  pep_df_mean = {
      "dRT": tdf['dRT'].mean(),
      "dPPM": tdf['dPPM'].mean(),
      "# Chargestates": tdf.groupby(['Date','Name'])['Chargestate'].nunique().mean(),
      "# Identified QC2 Peptides": tdf.groupby(['Date','Name'])['Peptide'].nunique().mean(),
  }
  setattr(df_widget, 'peps_mean', pep_df_mean)

def plot_calendar_hist(df_widget_value, df_widget:pn.widgets.Tabulator=df_widget):
  if (df_widget is None) or (type(df_widget.value != pd.DataFrame)) or (df_widget_df.shape[0]):
    return pd.DataFrame(columns=["Date","Value"]).hvplot.line(
          title="QC2 Run Yearly Distribution", 
          ylabel="# Runs",
          xlabel="Month of {}".format("N/A"))

  df_widget_df = df_widget.value
  calendar_hist = df_widget_df.assign(month=df_widget_df.Date.dt.month)\
        .groupby(by=['month']).count().hvplot.bar(y="Date")\
        .opts(title="QC2 Run Yearly Distribution", ylabel="# Runs",
              xlabel="Month of {}".format(next(iter(df_widget_df.Date.dt.year))))
  return calendar_hist

def update_ds(ds_key, mzqc_paths=mzqc_paths, df_widget:pn.widgets.Tabulator=df_widget, ds_descriptor:pn.pane.Markdown=ds_descriptor,
              date_range_slider:pn.widgets.DateRangeSlider=date_range_slider, checkbox_group:pn.widgets.CheckBoxGroup=checkbox_group):
  # sanity check
  if ds_key not in mzqc_paths.keys():
    print("ds_key not in mzqc_paths.")
    return

  if df_widget is None:
    print(">>> skipping enpty init", ds_key, df_widget)
    return

  print(">>> updating with {}".format(ds_key))

  # load new dataset
  main, tic, pep = load_ds(ds_key, mzqc_paths)
  # update df_widget

  df_widget.selection = list()
  setattr(df_widget, 'tics', tic)
  setattr(df_widget, 'peps', pep)
  df_widget.value = main
  update_ds_pep_stats(df_widget) 

  print(">>> updated df_widget")

  # update checkbox_group
  dir(checkbox_group)
  setattr(checkbox_group, 'options', main.columns.drop(['Date', 'RT range', 'MZ range', 'mzrange', 'Name']).to_list())  # do not deactivate unless all mzQC produce the same column layout
  setattr(checkbox_group, 'value', ['# MS1', '# MS2'])

  print(">>> updated chekcbox_group")

  # update date_range_slider
  date_range_slider.start=dt.datetime.combine(main.Date.min(),dt.datetime.min.time())
  date_range_slider.end=dt.datetime.combine(main.Date.max(),dt.datetime.min.time())
  date_range_slider.value=(date_range_slider.start,date_range_slider.end)

  print(">>> updated date_range_slider")

  ds_descriptor.value = descr_templ.format(s=str(main.Date.min().date()),
                                          e=str(main.Date.max().date()),
                                          p=dataset_mid_path[ds_key],
                                          n=str(len(main)),
                                          i="Thermo {}".format(next(iter(ds_key.split('_')))))

  print(">>> updated ds_descriptor")

  return


df_widget_max_select = 6
df_widget = pn.widgets.Tabulator(value=pd.DataFrame(columns=
  ['Date','# MS1','# MS2','# ID MS2','# Features','# ID Features', '# Signal fluct. ↓', '# Signal fluct. ↑','Name','mzrange',]),
                                  selectable='checkbox',
                                  show_index=False,
                                  pagination='local', page_size=10,
                                  disabled=True,
                                  hidden_columns=["Name",'mzrange'])
df_widget.selectable=df_widget_max_select
setattr(df_widget, 'tics', None)  # squirrel all ds_ into the df_widget? 
setattr(df_widget, 'peps', None)

ds_select = pn.widgets.Select(name='Select Dataset', options=list(mzqc_paths.keys()))

tics_pane = pn.bind(plot_tics, selection=df_widget.param.selection)

runsticker_pane = pn.bind(plot_runsticker, selection=df_widget.param.selection)

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

metrics_pane = pn.bind(plot_metrics_daterange, 
                       start=date_range_slider.param.value_start, end=date_range_slider.param.value_end,
                       selection=checkbox_group.param.value,)

cchart_pane = pn.bind(plot_ccharts, 
                       start=date_range_slider.param.value_start, 
                       end=date_range_slider.param.value_end,
                       )

descr_templ = """
# The Data Set

This dataset ({p}) is comprised of mzQC files from {n} runs,

which were measured over the timecourse between {s} and {e}, 

from QC2 samples on an instrument of type {i}.

You can in inspect the details in the following interactive charts.
"""
ds_descriptor = pn.pane.Markdown(descr_templ.format(s=str(dt.datetime.now().date()),
                                                               e=str(dt.datetime.now().date()),
                                                               p="PXD------",
                                                               n=str(len([])),
                                                               i="<Make> <Model>"))

calendar_hist_pane = pn.bind(plot_calendar_hist, df_widget_value=df_widget.param.value,)

ds_switch = pn.bind(update_ds, ds_key=ds_select.param.value)

row0 = pn.Row(ds_select)
row1 = pn.Row(ds_descriptor, calendar_hist_pane)
row2 = pn.Row(df_widget)
row3 = pn.Row(tics_pane, runsticker_pane)
internal_col = pn.Column(date_range_slider, checkbox_group)
row4 = pn.Row(internal_col, metrics_pane)
# row5 = pn.Row(tabs = pn.Tabs(cchart_pane, tabs_location='left'))   #,('Rbar control chart', cchart_pane)
col = pn.Column(row0, row1, row2, row3, row4, cchart_pane, ds_switch)  # bind must be included to be active
col.servable()



# ! python /home/vscode/.local/bin/panel serve crg_qc2_timeseries.py
# ! panel serve qc2_timeseries.py
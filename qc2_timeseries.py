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
pn.extension('tabulator')

# mzqc_basepath = "/content/drive/Shareddrives/mzqclib-manuscript/test data/CRG"
# mzqc_basepath = "mzqcs"
mzqc_basepath = "smalldevset"
dataset_mid_path = {"Lumos_2017": "PXD019888",
"Velos_2018": "PXD019889",
"Lumos_2018": "PXD019891",
"Velos_2017": "PXD019892"}
metadata_paths = {k: os.path.join(os.path.join(mzqc_basepath, v), "metadata_"+k+".json") for k,v in dataset_mid_path.items()}
mzqc_paths = {k: [os.path.join(os.path.join(mzqc_basepath, v), x) for x in os.listdir(os.path.join(mzqc_basepath, v)) if x.endswith(".mzqc")] for k,v in dataset_mid_path.items()}

ds_mzqcs = dict()
for dskey, dspaths in mzqc_paths.items():
  dskeyruns = dict()
  for mzqc_path in dspaths:
    with open(mzqc_path, "r") as file:
      mzqcobj = qc.JsonSerialisable.FromJson(file)
      dskeyruns[os.path.basename(mzqc_path)] = mzqcobj
  ds_mzqcs[dskey] = dskeyruns

def load_df(mzqc_dict: Dict[str, qc.MzQcFile]):
  df = pd.DataFrame([{'Date': next(iter([pd.to_datetime(cd.value) for cd in rowrunmzqc.runQualities[0].metadata.inputFiles[0].fileProperties if cd.accession=="MS:1000747"])),
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
                      'Name':rowrunname,
                      } for rowrunname, rowrunmzqc in mzqc_dict.items()])

  df.Date = df.Date.dt.date
  df.sort_values(by='Date', inplace = True) 
  return df

ds_dfs = {dskey: load_df(dskeyruns) for dskey, dskeyruns in ds_mzqcs.items()}
ds_select = pn.widgets.Select(name='Select Dataset', options=list(ds_dfs.keys()))
ds_select.value = next(iter(ds_dfs.keys()))  # this is important

df_widget_max_select = 6
df_widget = pn.widgets.Tabulator(value=next(iter(ds_dfs.values())), selectable='checkbox',
                                  show_index=False,
                                  pagination='local', page_size=10,
                                  disabled=True,
                                  hidden_columns=["Name",'mzrange'])
df_widget.selectable=df_widget_max_select
# TODO Does df_widget need to be servable?!?

def plot_tics(selection=[], ds_select=ds_select, ds_mzqcs=ds_mzqcs, ds_dfs=ds_dfs):
  if len(selection)==0:
    return pd.DataFrame(columns=["Retentiontime [s]","Relative Intensity of Ion Current"]).hvplot.line(title="Total Ion Chromatogram",
                                                            line_width=0.5, 
                                                            xlabel="Retentiontime [s]", 
                                                            ylabel="Relative Intensity of Ion Current",
                                                            frame_height=500,
                                                            frame_width=800,)
  
  data = ds_dfs[ds_select.value]
  with pn.param.set_values(col[1][1], loading=True):
    selection_names = data.iloc[selection].Name
    selection_dates = data.iloc[selection].Date

    selection_tics = [pd.DataFrame(next(iter([m.value for 
                                                m in ds_mzqcs[ds_select.value].get(n).runQualities[0].qualityMetrics if m.accession=="MS:4000104"]))) for 
                                                  n in selection_names]
    selection_tics = [tic.assign(**{"Date":str(d), "Name":n}) for n,d,tic in zip(selection_names,selection_dates,selection_tics)]

    tics_df =  pd.pivot(pd.concat(selection_tics), index="MS:1000894", columns="Date", values="MS:1000285")
    return tics_df.interpolate(method='linear').hvplot.line(legend='top_right', 
                                                            title="Total Ion Chromatogram", 
                                                            line_width=0.5, 
                                                            xlabel="Retentiontime [s]", 
                                                            ylabel="Relative Intensity of Ion Current",
                                                            frame_height=500,
                                                            frame_width=800,)

def plot_runsticker(selection=[], ds_select=ds_select, ds_dfs=ds_dfs):
  if len(selection)==0:
    selection_df = pd.DataFrame(columns=['Date','# MS1','# MS2','# ID MS2','# Features','# ID Features', '# Signal fluct. ↓', '# Signal fluct. ↑',])
    idax = selection_df[['# MS1','# MS2', '# ID MS2', 'Date']].set_index('Date').hvplot.barh(frame_height=200, frame_width=200)
    qaax = selection_df[['# Features','# ID Features', 'Date']].set_index('Date').hvplot.barh(frame_height=200, frame_width=200)
    mzax = selection_df.hvplot.line(x="Date", xlabel="m/z Range Setting", ylim=(100,1600),
                             invert=True, frame_height=200, frame_width=200).opts(yaxis='bare')
    flax = selection_df[['# Signal fluct. ↓', '# Signal fluct. ↑', 'Date']]\
            .set_index('Date').astype('int').hvplot.bar(rot=45, frame_height=200, frame_width=200)
    return (idax + qaax + mzax + flax).cols(1).opts(shared_axes=False)

  selection_df = ds_dfs[ds_select.value].iloc[selection]
  with pn.param.set_values(col[1][2], loading=True):
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

    return (idax + qaax + mzax + flax).cols(1).opts(shared_axes=False,)

tics_pane = pn.bind(plot_tics, selection=df_widget.param.selection)
runsticker_pane = pn.bind(plot_runsticker, selection=df_widget.param.selection)

date_range_slider = pn.widgets.DateRangeSlider(
    name='Date Range Slider',
    start=next(iter(ds_dfs.values())).Date.min(), end=next(iter(ds_dfs.values())).Date.max(),
    # value=(next(iter(ds_dfs.values())).Date.min()+dt.timedelta(days=10), 
    #       next(iter(ds_dfs.values())).Date.max()-dt.timedelta(days=10)),
    value=(next(iter(ds_dfs.values())).Date.min(), 
          next(iter(ds_dfs.values())).Date.max()),
    # step=24*3600*2*1000,
    # callback_throttle=100,
    tooltips=True,
)

checkbox_group = pn.widgets.CheckBoxGroup(
    name='Checkbox Group', value=['# MS1', '# MS2'], 
    options=next(iter(ds_dfs.values())).columns.drop('Date').to_list(),
    inline=False
)

def plot_metrics_daterange(start, end, selection, ds_select=ds_select, ds_dfs=ds_dfs):
  if not start or not end or not selection:
    return pd.DataFrame(columns=['Date','# MS1','# MS2',]).hvplot.line(x="Date", title="Date Range Metric Values", ylabel="Absolute Metric Values")
  try:
    with pn.param.set_values(col[3][2], loading=True):
      selection_df = ds_dfs[ds_select.value][selection+['Date']]
      selection_df = selection_df[ds_dfs[ds_select.value].Date.between(start.date(),end.date())]
      return selection_df.hvplot.line(x="Date", title="Date Range Metric Values", ylabel="Absolute Metric Values", )
  except Exception as e:
    print( "!!!exception",e )
    return pd.DataFrame(columns=['Date','# MS1','# MS2',]).hvplot.line(x="Date", title="Date Range Metric Values", ylabel="Absolute Metric Values")

metrics_pane = pn.bind(plot_metrics_daterange, 
                       start=date_range_slider.param.value_start, end=date_range_slider.param.value_end,
                       selection=checkbox_group.param.value,)

def update_ds(dskey, df_widget=df_widget, date_range_slider=date_range_slider, checkbox_group=checkbox_group):
  # update df_widget
  df_widget.value = ds_dfs[dskey]
  df_widget.selection = list()
  # update checkbox_group
  checkbox_group.options=ds_dfs[dskey].columns.drop('Date').to_list()
  checkbox_group.value=['# MS1', '# MS2']
  # update date_range_slider
  date_range_slider.start=ds_dfs[dskey].Date.min()
  date_range_slider.end=ds_dfs[dskey].Date.max()
  date_range_slider.value=(ds_dfs[dskey].Date.min(), ds_dfs[dskey].Date.max())

ds_switch = pn.bind(update_ds, dskey=ds_select.param.value)

row0 = pn.Row('# Column Explore', ds_select)
row1 = pn.Row('## Row1', tics_pane, runsticker_pane)
row2 = pn.Row('## Row2', df_widget)
internal_col = pn.Column(date_range_slider, checkbox_group)
row3 = pn.Row('## Row3', internal_col, metrics_pane)
col = pn.Column(row0, row1, row2, row3, ds_switch)  # bind must be included to be active
col.servable()



# !python /home/vscode/.local/bin/panel serve crg_qc2_timeseries.py
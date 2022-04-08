import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import pandas as pd
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from itertools import combinations
from Bio.Blast import NCBIXML
import glob
import plotly.graph_objects as go
import io
import numpy as np


def drawOrf(orfs, z=0, h=0.2):
    firstorf = np.min([x.location.start for x in orfs])
    trace_list = []
    for orf in orfs:
        if orf.type == 'CDS':
            start, stop, strand = orf.location.start, orf.location.end, orf.location.strand
            color = 'gray'
            linecolor = 'black'
            opacity=1
            x, y = draw_shape(start, stop, strand, z, h, firstorf)
            trace = go.Scatter(
                    x=x,
                    y=y,
                    marker=dict(size=3),
                    opacity=opacity,
                    fill='toself',
                    fillcolor=color,
                    line_color=linecolor,
                    line_width=3,
                    text=orf.qualifiers.get('product', 'na'),
                    hoverinfo='text')
            trace_list.append(trace)
    print(trace_list)
    return trace_list

def load_trna_trace(self,z=0,h=0.2):
    self.firstorf = np.min([x.location.start for x in self.orfs])
    self.lastorf = np.max([x.location.start for x in self.orfs])
    trace_list = []
    if len(self.tRNAs) > 0:
        for trna in self.tRNAs:
            print(trna)
            print(dir(trna.location))
            if trna.location is None:
                continue
            start, stop, strand = trna.location.start, trna.location.end, trna.location.strand
            linecolor = 'gold'
            opacity = 1
            x, y = draw_trna(start,stop,strand,z,h,self.firstorf,self.lastorf)
            trace = go.Scatter(x=x, y=y, marker=dict(size=1), opacity=opacity, fill='toself', fillcolor='gold', line_color='gold', text='{}|{}'.format(''.join(trna.product), self.annotations['organism']),hoverinfo='text')
            trace_list.append(trace)
    return trace_list

def load_repeat_trace(self,z=0,h=0.2):
        self.firstorf = np.min([x.location.start for x in self.orfs])
        self.lastorf = np.max([x.location.start for x in self.orfs])
        trace_list = []
        if len(self.repeat_regions) > 0:
            for rp in self.repeat_regions:
                start, stop, strand = rp.location.start, rp.location.end, rp.location.strand
                linecolor = 'pink'
                opacity = 1
                x, y = draw_repeat(start,stop,strand,z,h,self.firstorf,self.lastorf)
                trace = go.Scatter(x=x, y=y, marker=dict(size=1), opacity=opacity, fill='toself', fillcolor='pink', line_color='pink', text='{}|{}'.format(''.join(rp.rpt_family), self.annotations['organism']),hoverinfo='text')
                trace_list.append(trace)
        return trace_list
      

def draw_shape(start, stop, strand, z, h, firstorf):
    start, stop = int(start), int(stop)
    arrow_offset = 150
    if strand == 1:
        x=(start, start+arrow_offset, start, stop-arrow_offset, stop, stop-arrow_offset, start)
        y=(z, z+h/2, z+h, z+h, z+h/2, z, z)
    elif strand == -1:
        x=(start+arrow_offset, start, start+arrow_offset, stop, stop-arrow_offset, stop, start+arrow_offset)
        y=(z, z-h/2, z-h, z-h, z-h/2, z, z)
    else:
        z = z+0.5 #offset height
        x=(start, start, stop, stop, start)
        y=(z, z+h, z+h, z, z)
    return x-firstorf,y

def draw_trna(start,stop,strand,z,h,firstorf,lastorf):
    start,stop = int(start), int(stop)
    if strand == 1:
        x=(start, start, stop, stop, start)
        y=(z, z+h, z+h, z, z)
    else:
        x=(start, start, stop, stop, start)
        y=(z, z-h, z-h, z, z)
    return x,y

def draw_repeat(start,stop,strand,z,h,firstorf,lastorf):
    start,stop = int(start), int(stop)
    if strand == 1:
        x=(start,start,stop,stop,start)
        y=(z,z+h,z+h,z,z)
    else:
        x=(start,start,stop,stop,start)
        y=(z,z-h,z-h,z,z)
    return x,y

def graphing(traces):
    fig = go.Figure(layout={'width':1200,'height':1200})
    fig.update_yaxes(range=[-5, 5])
    [fig.add_trace(x) for x in traces] 
    fig.layout.plot_bgcolor = 'white'
    fig.layout.paper_bgcolor = 'white'
    fig.update_layout(showlegend=False)
    fig.write_html('/hdd-roo/DHS/figures/file.html')
    


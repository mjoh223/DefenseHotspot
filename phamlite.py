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

def colorofsystem(sys):
    colors = [
        '#993c99',
        '#a53991',
        '#af3788',
        '#b73780',
        '#be3877',
        '#c33a6e',
        '#c73e65',
        '#c9435d',
        '#ca4a55',
        '#cb504d',
        '#ca5746',
        '#c85e40',
        '#c5663a',
        '#c26d35',
        '#bd7432',
        '#b87b2f',
        '#b3812e',
        '#ad872f',
        '#a78e32',
        '#a09336',
        '#99993c']
    systems = [
        'RM',
        'Shedu',
        'BREX',
        'Gabija',
        'Zorya',
        'SspBCDE',
        'DarTG',
        'Septu',
        'Gao_Qat',
        'Druantia',
        'Dnd',
        'AbiEii',
        'Lamassu',
        'Gao_RL',
        'Wadjet',
        'Hachiman',
        'DISARM',
        'CBASS',
        'Gao_Upx',
        'Kiwa',
        'Thoeris',]
    mydict = dict(zip(systems,colors))
    return mydict[sys]

def drawOrf(systems, contig, assembly, orfs, z=0, h=0.2):
    try:
        firstorf = np.min([x.location.start for x in orfs])
    except:
        firstorf = 0
    trace_list = []
    sys_list = []
    for i, orf in enumerate(orfs):
        #if i < 7:
        #    color = '#3c9999'
        #    opacity = 1
        #elif i > len(orfs)-7:
        #    color = '#3c9999'
        #    opacity = 1
        #else:
        #    color = 'gray'
        #    opacity = 0.4
        if orf.type == 'CDS' and 'protein_id' in orf.qualifiers:
            gembase_id = '{}q{}_{}'.format(assembly,contig,orf.qualifiers['protein_id'][0].split('.')[0].replace('_',''))
            print(gembase_id)
            #wpid = orf.qualifiers['protein_id'][0].replace('_','').split('.')[0]
            start, stop, strand = orf.location.start, orf.location.end, orf.location.strand
            color = 'gray'
            opacity = 0.4
            linecolor = 'black'
            x, y = draw_shape(start, stop, strand, z, h, firstorf)
            sys = ''
            if gembase_id in systems:
                sys = systems[gembase_id]
                print(sys)
                sys_list.append(sys)
                color = colorofsystem(sys)
                print(color)
                #color = '#3c9999'
                opacity = 1
            product = orf.qualifiers.get('product', ['unknown'])[0]+'_'+sys
            trace = go.Scatter(
                    x=x,
                    y=y,
                    marker=dict(size=2),
                    opacity=opacity,
                    fill='toself',
                    fillcolor=color,
                    line_color=linecolor,
                    line_width=2,
                    text=gembase_id,
                    hoverinfo='text')
            trace_list.append(trace)
        
        if orf.type == 'tRNA':
            print(orf.location)
            start = orf.location.start
            stop = orf.location.end
            strand = orf.location.strand
            x, y = draw_trna(start, stop, strand, z, h, firstorf )
            color = 'purple'
            trace = go.Scatter(
                    x=x,
                    y=y,
                    marker=dict(size=1),
                    opacity=1,
                    fill='toself',
                    fillcolor=color,
                    line_color='purple',
                    line_width=1,
                    text='tRNA',
                    hoverinfo='text')
            trace_list.append(trace)
    sys_str = '_'.join(sys_list)
    sys_str = sys_str.replace(' ', '_')
    return trace_list, sys_str

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
        y=(h/2+z, h/2+z-h/2, h/2+z-h, h/2+z-h, h/2+z-h/2, h/2+z, h/2+z)
    else:
        z = z+0.5 #offset height
        x=(start, start, stop, stop, start)
        y=(z, z+h, z+h, z, z)
    return x-firstorf, y

def draw_trna(start,stop,strand,z,h,firstorf):
    start,stop = int(start), int(stop)
    if strand == 1:
        x=(start, start, stop, stop, start)
        y=(z, z+h, z+h, z, z)
    else:
        x=(start, start, stop, stop, start)
        y=(h/2+z, h/2+z-h, h/2+z-h, h/2+z, h/2+z)
    return x-firstorf, y

def draw_repeat(start,stop,strand,z,h,firstorf,lastorf):
    start,stop = int(start), int(stop)
    if strand == 1:
        x=(start,start,stop,stop,start)
        y=(z,z+h,z+h,z,z)
    else:
        x=(start,start,stop,stop,start)
        y=(z,z-h,z-h,z,z)
    return x,y

def graphing(traces, contig, size, sys):
    size = np.absolute(size)
    fig = go.Figure(layout={'width':size/50,'height':400})
    fig.update_yaxes(range=[-3, 3], visible=False)
    fig.update_xaxes(range=[-1000, size+1000]) #was 45000
    [fig.add_trace(x) for x in traces] 
    fig.layout.plot_bgcolor = 'white'
    fig.layout.paper_bgcolor = 'white'
    fig.update_layout(showlegend=False, margin=dict(l=0,r=0,b=0,t=0))
    fig.write_image('/hdd-roo/DHS/DHS_df_svg/{}_{}_{}.svg'.format(size, contig, sys))
    fig.write_html('/hdd-roo/DHS/DHS_df_html/{}_{}_{}.html'.format(size, contig, sys))
    


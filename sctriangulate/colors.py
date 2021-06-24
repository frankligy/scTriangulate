from matplotlib import cm
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex
import copy

'''
background greyed colormap 
'''
def bg_greyed_cmap(cmap_str):
    # give a matplotlib cmap str, for instance, 'viridis' or 'YlOrRd'
    cmap = copy.copy(cm.get_cmap(cmap_str))
    cmap.set_under('lightgrey')
    return cmap



# zeileis_28 was taken from scanpy: https://github.com/theislab/scanpy/blob/master/scanpy/plotting/palettes.py
# and they noted the original source as below:
# https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# update 1
# orig reference http://epub.wu.ac.at/1692/1/document.pdf

zeileis_28 = [
    "#023fa5",
    "#7d87b9",
    "#bec1d4",
    "#d6bcc0",
    "#bb7784",
    "#8e063b",
    "#4a6fe3",
    "#8595e1",
    "#b5bbe3",
    "#e6afb9",
    "#e07b91",
    "#d33f6a",
    "#11c638",
    "#8dd593",
    "#c6dec7",
    "#ead3c6",
    "#f0b98d",
    "#ef9708",
    "#0fcfc0",
    "#9cded6",
    "#d5eae7",
    "#f3e1eb",
    "#f6c4e1",
    "#f79cd4",
    # these last ones were added:
    '#7f7f7f',
    "#c7c7c7",
    "#1CE6FF",
    "#336600",
]

# godsnot_102 was taken from scanpy: https://github.com/theislab/scanpy/blob/master/scanpy/plotting/palettes.py
# the author noted the original source as below:
# from http://godsnotwheregodsnot.blogspot.de/2012/09/color-distribution-methodology.html
godsnot_102 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00",
    "#1CE6FF",
    "#FF34FF",
    "#FF4A46",
    "#008941",
    "#006FA6",
    "#A30059",
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
    "#5A0007",
    "#809693",
    "#6A3A4C",
    "#1B4400",
    "#4FC601",
    "#3B5DFF",
    "#4A3B53",
    "#FF2F80",
    "#61615A",
    "#BA0900",
    "#6B7900",
    "#00C2A0",
    "#FFAA92",
    "#FF90C9",
    "#B903AA",
    "#D16100",
    "#DDEFFF",
    "#000035",
    "#7B4F4B",
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
    "#404E55",
    "#0089A3",
    "#CB7E98",
    "#A4E804",
    "#324E72",
]

pub_icgs2 = [
    '#F26D6D',  # red
    '#BF9004',  # brown
    '#62BF04',  # blue
    '#2BB7EC',  # cyan
    '#A38BFD',  # purple
    '#F263DA',  # pink
]

pub_shap = [
    '#F2075D',   # red
    '#158BFB',    # blue
]


def pick_n_colors(n):
    if n <= 10:
        _colors = [to_hex(color) for color in cm.get_cmap('tab10').colors[:n]]
    elif n > 10 and n <= 20:
        _colors = [to_hex(color) for color in cm.get_cmap('tab20').colors[:n]]
    elif n > 20 and n <= 28:
        _colors = zeileis_28[:n]
    elif n > 28 and n <= 102:
        _colors = godsnot_102[:n]
    elif n > 102:
        _colors = [colors.to_hex(cm.jet(round(i))) for i in np.linspace(0,255,n)]
    return _colors

def colors_for_set(setlist):  # a list without redundancy
    length = len(setlist)
    _colors = pick_n_colors(n=length)
    cmap = pd.Series(index=setlist,data=_colors).to_dict()
    return cmap

'''
below stores the nice cmap I encoutered in my research
'''

# Nathan's Yellow-blue schema
# yellow-blue colormap
cdict = {
    'red':((0.0,0.0,0.0),
           (0.5,0.0,0.0),
           (1.0,1.0,1.0)),
    'green':((0.0,0.8,0.8),
             (0.5,0.0,0.0),
             (1.0,1.0,1.0)),
    'blue':((0.0,1.0,1.0),
            (0.5,0.0,0.0),
            (1.0,0.0,0.0))
}

ywb_cmap = LinearSegmentedColormap('yellow_blue',segmentdata=cdict)



# SHAP pink-blue schema
cdict = {'red':((0.0,0.0,0.0),
                (1.0,1.0,1.0)),
         'green':((0.0,0.5,0.5),
                  (0.73,0.0,0.0),
                  (1.0,0.0,0.0)),
         'blue':((0.0,1.0,1.0),
                 (1.0,0.0,0.0))}
pwb_cmap = LinearSegmentedColormap('shap', segmentdata=cdict)


# scPhere confusion matrix schema
cdict = {'red':((0.0,0.43,0.43),   # red chrome is 0.43 around (both left and right) 0.0 point, then increase to
                (0.45,1.0,1.0),    # 1.0 around 0.45 point, finally arrive
                (1.0,0.95,0.95)),  # 0.95 around 1.0 point
         'green':((0.0,0.61,0.61),
                  (0.45,1.0,1.0),
                  (1.0,0.27,0.27)),
          'blue':((0.0,0.85,0.85),
                  (0.4,0.96,0.96),
                  (1.0,0.18,0.18))}
scphere_cmap = LinearSegmentedColormap('scphere', segmentdata=cdict)







    

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.collections import LineCollection
import matplotlib as mpl
import numpy as np
import pylab
import os
import pickle

def get_data(filename, cut = 0):
    
    print("reading file...")

    f1 = ""
    f2 = filename
    f3 = "_QEP010031_1.txt"
    seq = f1 + f2 + f3
    msec = []
    data = []
    file = open(seq)
    row  = 0
    for i in file:
	if row/1000.0 == int(row/1000):
	    print str(row/1000)
        if row == 0 and "Begin Spectral Data" in i:
            row = 1
        elif row == 1:
            nm = [float(j) for j in i.split()]
            row += 1
        elif row > 1:
            wd = i.split()
            ms = int(wd[1])
            ab = [float(wd[j]) for j in range(2, len(wd))]
            assert len(nm) == len(ab)
            msec.append(ms)
            data.append(ab)
            row += 1
            
        if row != 0 and row > cut:
            break
    file.close()
    
    return (data, nm)

def get_wavelength(nm, wl):
    """returns the index of the wavelength provided"""
    ret = 0
    while wl > nm[ret]:
        ret += 1
        
    return ret
    
def get_column(data, wl):
    ret = []
    
    for i in data:
        ret.append(i[wl])
        
    return ret

def to2D(single, dim, wl):
    x = dim[0]
    y = dim[1]    
    if(x*y != len(single)):
        print("dimension error")
        return None
    r = -1

    ret = [[0 for i in range(x)] for j in range(y)]
    
    for i in range(0,len(single)):
        j = i % x
        if j == 0:
            r+=1
        ret[r][j] = single[i][wl]
        
    return ret

def to2D2(single, x, wl):
    r = -1

    ret = [[] for i in range(int(len(single)/x)+1)]
    
    for i in range(0,len(single)):
        j = i % x
        if j == 0:
            r+=1
        print("r: %d, j: %d, i: %d" % (r,j,i))
        ret[r].append(single[i][wl])
        
    while len(ret[r]) < len(ret[0]):
        ret[r].append(0)
        
    print("x: %d, y: %d" % (len(ret[0]), len(ret)))
    return ret
    
def show_image(image_array, colorbar=True, amp=False, vmin=None, vmax=None, show=True, name=None, dir=None):
    fig = plt.figure()
    ax = plt.gca()
    
    if(amp):
        cmap = "RdBu"
        vmax = np.amax(image_array)
        vmin = -vmax
        img = ax.imshow(image_array, cmap=cmap,vmin=vmin,vmax=vmax)
    else:
        cmap = "Blues"
        vmax = np.amax(image_array)
        vmin = 0
        img = ax.imshow(image_array, cmap=cmap)
        
    if(colorbar):
        cbar = fig.colorbar(img, ticks=[vmin, (vmax+vmin)/2, vmax])
    if(name != None):
        save(name, dir)
    if(show):
        pylab.show()
    
def scatter_plot(data, xlog=False, ylog=False, highlight=0, default=True, left=None, right=None, bottom=None, top=None, show=True, name=None, dir=None):
    x = data[0]
    y = data[1]
    fig = plt.figure()
    ax = plt.gca()
    if not default:
        if(left == None):
            left = min(x)
        if(right == None):
            right = max(x)
        if(bottom == None):
            bottom = min(y)
        if(top == None):
            top = max(y)
        ax.set_xlim(left=left,right=right)
        ax.set_ylim(bottom=bottom,top=top)
        
    if highlight == 0:
        ax.scatter(x, y, marker=".")
    else:
        ax.scatter(x[:highlight], y[:highlight], marker="o")
        ax.scatter(x[highlight:], y[highlight:], marker=".")
    if(xlog):
        ax.set_xscale('log')
    if(ylog):
        ax.set_yscale('log')
    if(name != None):
        save(name, dir)
    if(show):
        pylab.show()
        
def save(name, dir=None):       
    if(dir == None):
        os.chdir("./figures")
        plt.savefig(name)
        os.chdir("..")
    else:
        os.chdir("./figures")
        try:
            os.mkdir("./" + dir)
        except:
            pass
        os.chdir("./" + dir)
        plt.savefig(name)
        os.chdir("../..")
    
def average_spectrum(data, nm):
    avg = [0 for i in nm]
    
    for i in range(0,len(nm)):
        temp = []
        for j in data:
            temp.append(j[i])
        avg[i] = np.mean(temp)    
        
    return avg

def absorbance(data, nm):
    avg = average_spectrum(data, nm)
    absorbance = [[0 for x in nm] for y in data]
    for i in range(len(data)):
        for j in range(len(nm)):
            absorbance[i][j] = 0
            #print(data[i][j])
            if(data[i][j] <= 0 or avg[j] <= 0):
                continue
            #print(i, "/", len(data))
            absorbance[i][j] = np.log10(data[i][j]/avg[j])
            
    return absorbance

def extract_components(matrix, comp, column=True):
    if column:
        ret = [i[:comp] for i in matrix]
        return ret
    else:
        ret = matrix[:comp]
        return ret

def reverse_svd(U, W, Vt, comp):    
    newU = extract_components(U, comp, column=True)
    
    newVt = extract_components(Vt, comp, column=False)
    
    newW = [[0 for x in range(comp)] for y in range(comp)]
    
    for i in range(comp):
        newW[i][i] = W[i]
    
    ret = np.dot(np.dot(newU, newW), newVt)
    
    return ret

def versus_graph(matrix, t, indexes, cmap='gist_rainbow', column=True, xlog=False, ylog=False, left=None, right=None, bottom=None, top=None, margin=None, gradient=False, show=True, name=None, dir=None):
    if column:
        x = [i[indexes[0]] for i in matrix]
        y = [i[indexes[1]] for i in matrix]
    else:
        x = matrix[indexes[0]]
        y = matrix[indexes[1]]
    
    fig = plt.figure()
    ax = plt.gca()
    plt.ylabel(str(indexes[1]))
    plt.xlabel(str(indexes[0]))
    
    if gradient:
        lc = gradient_line((x,y), t)
        ax.add_collection(lc)
    else:
        plt.plot(x, y)
        plt.plot([x[0], x[-1]], [y[0], y[-1]], 'ro')
    
    if(left == None):
        left = min(x)
    if(right == None):
        right = max(x)
    if(bottom == None):
        bottom = min(y)
    if(top == None):
        top = max(y)
    
    if margin != None:
        ax.set_xlim(left=(left-margin),right=(right+margin))
        ax.set_ylim(bottom=(bottom-margin),top=(top+margin))
    else:
        ax.set_xlim(left=left, right=right)
        ax.set_ylim(bottom=bottom, top=top)
    
    if show:
        pylab.show()
        
def gradient_line(funcs, t, cmap='gist_rainbow', width=3):
    x = funcs[0]
    y = funcs[1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    #points = np.array([x,y])
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    lc = LineCollection(segments, cmap=plt.get_cmap(cmap), norm=plt.Normalize(min(t), max(t)))
    lc.set_array(t)
    lc.set_linewidth(width)
    
    return lc

dark_threshold = 15

def is_dark(spec, threshold=dark_threshold):
    return np.mean(spec) <= threshold

def align_dark(data, threshold=dark_threshold, min=3):
    """starting with first dark spectrum"""
    print("aligning by dark spectrums...")
    i = 0
    
    while i < len(data) and not is_dark(data[i], threshold):
        i+=1
    
    c = 0
    while i < len(data) and is_dark(data[i], threshold):
        i+=1
        c+=1
        
    ret = []
    j = 0
    
    while i < len(data):
        ret.append([])   
        while i < len(data) and not is_dark(data[i], threshold):
            ret[j].append(data[i])
            i+=1
            
        c = 0
        while i < len(data) and is_dark(data[i], threshold):
            i+=1
            c+=1
        j += 1
            
    return ret

def isolate_wavelength(twoD_image, wl):
    print("isolating single wavelength...")
    ret = [[] for i in twoD_image]
    
    for i in range(len(twoD_image)):
        for j in twoD_image[i]:
            ret[i].append(j[wl])
            
    return ret

def square(twoD):
    print("squaring image...")
    
    twoD_image = twoD[:(len(twoD)-1)]
    
    ret = [[] for i in twoD_image]
    
    min = -1
    
    for i in twoD_image:
        if min == -1 or len(i) < min:
            min = len(i)
            
    for i in range(len(twoD_image)):
        ret[i] = twoD_image[i][:min]
        
    print(min)
    return ret

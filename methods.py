import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.collections import LineCollection
import matplotlib as mpl
import numpy as np
import pylab
import os
import pickle

def get_data(filename):
    """
    for root, dirs, files in os.walk("./pickle"):
      for file in files:
        if file.endswith(filename + ".data.p"):
          print(filename + ".data.p loading...")
          os.chdir("pickle")
          with open(filename + ".data.p", "rb") as f:
            v1 = pickle.load(f)
          os.chdir("..")
          return (v1[0], v1[1])
    """
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
    file.close()
    """    
    v1 = (data, nm)

    os.chdir("pickle")
    with open(filename + ".data.p", "wb") as f:
      pickle.dump(v1, f)
    os.chdir("..")
    """
    return (data, nm) #should be recursive call to self???
    
def get_column(data, wl):
    ret = []
    
    for i in data:
        ret.append(i[wl])
        
    return ret

def square(data, y, filename=None):

  if filename is not None:
    for root, dirs, files in os.walk("./pickle"):
      for file in files:
        if file.endswith(filename + ".square.p"):
          print(filename + ".square.p loading...")
          os.chdir("pickle")
          with open(filename + ".square.p", "rb") as f:
            v1 = pickle.load(f)
          os.chdir("..")
          return v1
  
  x = int(len(data)/y)

  #tolerance
  t = 10
  #cursor
  c = 0
  #break spots
  breaks = []

  dt = [np.mean(i) for i in data]

  for i in range(y):
    c += x

    #list of corrcoef
    cc = []
    best = [-1,-1]

    #meme code???
    for j in range(t):
      for k in [-1,1]:
        v1 = t*k
        v2 = c-x+v1
        if v2 < 0: v2 = 0

        t1 = dt[v2:c+v1]
        t2 = dt[c+v1:c+len(t1)+v1]

        t3 = np.mean(np.corrcoef(t1, t2)[0])

        print(t3)

        cc.append(t3)
        if best[1] == -1 or abs(1 - t3) < best[1]:
          best[1] = abs(1 - t3)
          best[0] = j
    
    breaks.append(best[0])

  temp = len(data) - 1
  data2D = [0 for i in range(y)]

  for i in reversed(range(len(breaks))):
    data2D[i] = data[breaks[i]:temp]

    temp = breaks[i]

  ret = [data2D[i][:x] for i in range(y)]

  if filename is not None:
    os.chdir("pickle")
    with open(filename + ".square.p", "wb") as f:
      pickle.dump(ret, f)
    os.chdir("..")

  return ret

def to2D(single, dim, wl):
    x = dim[0]
    y = dim[1]    
    if(x*y != len(single)):
        print("dimention error")
        return None
    r = -1

    ret = [[0 for i in range(x)] for j in range(y)]
    
    for i in range(0,len(single)):
        j = i % x
        if j == 0:
            r+=1
        ret[r][j] = single[i][wl]
        
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
    

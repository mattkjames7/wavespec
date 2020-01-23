# wavespec
Some spectral analysis tools for analyzing waves in data.

## Installation

### Using ```pip3```:

```pip3 install wavespec --user```

### Installing the wheel using ```pip3```:

```pip3 install wavespec-0.0.1-py3-none-any.whl --user```

### From git:
 
```
git clone https://github.com/mattkjames7/wavespec
cd wavespec
python3 setup.py install --user
```

## Usage

```python
import wavespec as ws
```

### Fast Fourier Transform (FFT)

```python
power,phase,freq,fr,fi = ws.Fourier.FFT(t,x,WindowFunction=None,Param=None)
```

### Lomb-Scargle (LS)

```python
P,A,phi,a,b = ws.LombScargle.LombScargle(t,x0,f,Backend='C++',WindowFunction=None,Param=None)
```

### Spectrograms

```python
Nw,LenW,Freq,out = ws.Spectrogram.Spectrogram(t,v,wind,slip,Freq=None,Method='FFT',WindowFunction=None,Param=None,Detrend=True,FindGaps=True,GoodData=None,Quiet=True,LenW=None)

ax,Nw,LenW,Freq,Spec = ws.Spectrogram.PlotSpectrogram(t,v,wind,slip,Freq=None,Method='FFT',WindowFunction=None,Param=None,Detrend=True,FindGaps=True,GoodData=None,Quiet=True,LenW=None,fig=None,maps=[1,1,0,0],PlotType='Pow',scale=None,zlog=False,TimeAxisUnits='s',FreqAxisUnits='Hz')
```

### 3D Spectrograms

```python
Nw,LenW,Freq,Spec = ws.Spectrogram.Spectrogram3D(t,vx,vy,vz,wind,slip,Freq=None,Method='FFT',WindowFunction=None,Param=None,Detrend=True,FindGaps=False,GoodData=None)
```

### Tests

```python
ws.Test.TestLS(A=[1.0,2.0],f=[0.04,0.1],phi=[0.0,90.0],Backend='C++')
ws.Test.TestPolarization(xPow=2.0,xPhase=0.0,yPow=1.0,yPhase=40.0)
ws.Test,TestSpectrogram()
```

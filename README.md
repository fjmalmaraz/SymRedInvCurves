# *SymRedInvCurves*:
# Software for continuation of Symmetric Invariant curves for Skew-Product 

This is a software to find 1--parameter families of symmetric invariant curves of a forced 1-D
maps (a.k.a skew-product)  with Lyapunov exponent constant. This software has been used to perform  the numerical continuation for the paper 
*On non-smooth pitchfork bifurcations in invertible quasi-periodically forced 1-D maps* by [A. Jorba](http://www.maia.ub.es/~angel/) (U. Barcelona), [F.J. Munoz-Almaraz](https://www.uchceu.es/directorio/fj-munoz-almaraz) (U. CEU-Cardenal Herrera) and J.C. Tatjer (U. Barcelona) . 

## Software dependencies 

- A C compiler, GCC-4.2 compiler has been used for all the calculations in the paper 

- Fast Fourier Transform West [FFTW](http://www.fftw.org/) with a public repository in Github [https://github.com/FFTW/fftw3](https://github.com/FFTW/fftw3).

## Run instructions 

Execute the following commands on a terminal 

```sh
$ git clone https://github.com/fjmalmaraz/SymRedInvCurves SymRedInvCurves
$ cd SymRedInvCurves
$ make
$ ./cont_complex.exe
```

## Settings 

The forced 1-D map can be change modifying the file [cont_complex.c](https://github.com/fjmalmaraz/SymRedInvCurves/blob/master/cont_complex.c). Parameters to control continuation method can be modified in [parameter.h](https://github.com/fjmalmaraz/SymRedInvCurves/blob/master/parameter.h). A description is provided in the [documetation](https://github.com/fjmalmaraz/SymRedInvCurves/blob/master/doc/SymRedInvCurves.pdf) (doc folder).

## Licence (BSD) 

> Copyright (c) 2016, Francisco Javier Muñoz-Almaraz
> All rights reserved.
>
> 
> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are met:
>
> 1. Redistributions of source code must retain the above copyright
> notice, this list of conditions and the following disclaimer.
>
> 2. Redistributions in binary form must reproduce the above copyright
> notice, this list of conditions and the following disclaimer in the
> documentation and/or other materials provided with the distribution.
>
> 3. All advertising materials mentioning features or use of this software
> must display the following acknowledgement:
> This product includes software developed by Francisco Javier Muñoz-Almaraz.
>
> 4. Neither the name of the Francisco Javier Muñoz-Almaraz nor the
> names of its contributors may be used to endorse or promote products
> derived from this software without specific prior written permission.
> 
> THIS SOFTWARE IS PROVIDED BY Francisco Javier Muñoz-Almaraz ''AS IS'' AND ANY
> EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
> WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
> DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
> DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
> (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
> LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
> ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
> (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
> SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

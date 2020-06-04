
# **The Landau-Lifshitz-Gilbert equation (LLG)**
<br>

## The main equation<br>
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;M_{s}" title="M_{s}" /> を飽和磁化、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\gamma&space;_{0}" title="\gamma _{0}" /> を磁気回転比、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\alpha" title="\alpha" /> を減衰定数とすると、
ドメイン構造は局所磁化
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{M(r)}" title="\mathbf{M(r)}" /> の空間分布によって記述されるため、次のLLG方程式を用いて表すことができる。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\large&space;(1&plus;\alpha^{2}&space;)\frac{\partial&space;\mathbf{M}}{\partial&space;t}&space;=&space;-\gamma&space;_{0}\mathbf{M}\times&space;\mathbf{H_{eff}}&space;-&space;\frac{\gamma_{0}\alpha}{M_{s}}&space;\mathbf{M}\times&space;\left&space;(&space;\mathbf{M}\times&space;\mathbf{H_{eff}}&space;\right&space;)\:&space;\:&space;\:&space;\:&space;\:&space;(1)"><br><br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{eff}}"> は有効磁界であり、次のように表すことができる。
（<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mu&space;_{0}" title="\mu _{0}"> は真空の透磁率）<br><br>

<img src="https://latex.codecogs.com/gif.latex?\mathbf{H_{eff}}=-\frac{1}{\mu_{0}}\frac{\partial&space;E}{\partial\mathbf{M}}=-\frac{1}{\mu_{0}}\frac{\left&space;(&space;E_{anis}&plus;E_{exch}&plus;E_{ms}&plus;E_{external}&plus;E_{elastic}&space;\right&space;)}{\partial\mathbf{M}}\:&space;\:&space;\:&space;\:&space;\:&space;(2)"><br><br>

このプログラム内では, 簡単のため 
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{M}=M_{s}\mathbf{m}"> とする。
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{anis}"> は結晶磁気異方性エネルギーであり、立方晶では<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{anis}=\int&space;\left&space;[&space;K_{1}&space;\left&space;(&space;m_{1}^{2}m_{2}^{2}&plus;m_{1}^{2}m_{3}^{2}&plus;m_{2}^{2}m_{3}^{2}&space;\right&space;)&plus;K_{2}m_{1}^{2}m_{2}^{2}m_{3}^{2}\right&space;]dV\:&space;\:&space;\:&space;\:&space;\:&space;(3)"><br><br>

となる。この時、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;m_{i}"> は単位磁化ベクトルの成分であり、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;K_{1},&space;K_{2}"> は異方性定数を表している。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{exch}"> は交換エネルギーであり、磁化の向きの空間変動によってのみ決定され、<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{exch}=A\int\left&space;(&space;m_{1,1}^{2}&plus;m_{1,2}^{2}&plus;m_{1,3}^{2}&plus;&space;m_{2,1}^{2}&plus;m_{2,2}^{2}&plus;m_{2,3}^{2}&space;&plus;&space;m_{3,1}^{2}&plus;m_{3,2}^{2}&plus;m_{3,3}^{2}&space;\right&space;)dV\:&space;\:&space;\:&space;\:&space;\:&space;(4)"><br><br>

となる。この時、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;A"> は単位磁化ベクトルの成分である。また、添え字内のコンマは空間微分を表しており、以降
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;m_{i,j}=\frac{\partial&space;m_{i}}{\partial&space;x_{j}}"> とし、
<img src="https://latex.codecogs.com/gif.latex?x_{j}"> をデカルト座標の j 番目の位置ベクトル成分とする。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{ms}"> は静磁気エネルギーであり、<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{ms}=-\frac{1}{2}\mu&space;_{0}M_{s}\int&space;\mathbf{H_{d}\cdot&space;m}dV\:&space;\:&space;\:&space;\:&space;\:&space;(5)" ><br><br>

となる。この時、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{d}}"> は磁気モーメント間の長距離相互作用によって決定される浮遊磁場であり、

<img src="https://latex.codecogs.com/gif.latex?H_{d1,1}&plus;H_{d2,2}&plus;H_{d3,3}=-M_{s}\left&space;(&space;m_{d1,1}&plus;m_{d2,2}&plus;m_{d3,3}&space;\right&space;)\:&space;\:&space;\:&space;\:&space;\:&space;(6)" ><br><br>

と表すことができる。ここで、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{di}}"> は
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{d}}"> の構成要素である（添え字内のコンマは空間微分）。また、<br><br>

<img src="https://latex.codecogs.com/gif.latex?H_{di}=-{\phi&space;}_j&space;\:&space;\:&space;\:&space;\:&space;\:&space;(7)"><br><br>

のように磁気スカラーポテンシャル
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{\phi_{j}}"> を導入することで、(6)式は<br><br>

<img src="https://latex.codecogs.com/gif.latex?\Delta&space;\phi&space;=M_s(m_{1,1}&plus;m_{2,2}&plus;m_{3,3})\:&space;\:&space;\:&space;\:&space;\:&space;(8)" >

のように書き直される。

[](
式番号の付け方\
<<\:&space;\:&space;\:&space;\:&space;\:&space;(7)>>
)
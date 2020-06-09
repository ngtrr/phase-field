
# **The Landau-Lifshitz-Gilbert equation (LLG)**
<br>

## 方程式の概要<br>
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;M_{s}" title="M_{s}" /> を飽和磁化、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\gamma&space;_{0}" title="\gamma _{0}" /> を磁気回転比、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\alpha" title="\alpha" /> をダンピング定数とすると、
ドメイン構造は局所磁化
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{M(r)}" title="\mathbf{M(r)}" /> の空間分布によって記述されるため、次のLLG方程式を用いて表すことができる。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\large&space;(1&plus;\alpha^{2}&space;)\frac{\partial&space;\mathbf{M}}{\partial&space;t}&space;=&space;-\gamma&space;_{0}\mathbf{M}\times&space;\mathbf{H_{eff}}&space;-&space;\frac{\gamma_{0}\alpha}{M_{s}}&space;\mathbf{M}\times&space;\left&space;(&space;\mathbf{M}\times&space;\mathbf{H_{eff}}&space;\right&space;)\:&space;\:&space;\:&space;\:&space;\:&space;(1)"><br><br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{eff}}"> は有効磁界であり、次のように表すことができる。
（<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mu&space;_{0}" title="\mu _{0}"> は真空の透磁率）<br><br>

<img src="https://latex.codecogs.com/gif.latex?\mathbf{H_{eff}}=-\frac{1}{\mu_{0}}\frac{\partial&space;E}{\partial\mathbf{M}}=-\frac{1}{\mu_{0}}\frac{\left&space;(&space;E_{anis}&plus;E_{exch}&plus;E_{ms}&plus;E_{external}&plus;E_{elastic}&space;\right&space;)}{\partial\mathbf{M}}\:&space;\:&space;\:&space;\:&space;\:&space;(2)"><br><br>

エネルギー項は順に結晶磁気異方性、交換、静磁、Zeeman、弾性エネルギーであり、下記に求め方を順に示す。また、このプログラム内では簡単のため、 
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{M}=M_{s}\mathbf{m}"> とする。<br><br>

## 結晶磁気異方性エネルギー<br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{anis}"> は結晶磁気異方性エネルギーであり、立方晶では<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{anis}=\int&space;\left&space;[&space;K_{1}&space;\left&space;(&space;m_{1}^{2}m_{2}^{2}&plus;m_{1}^{2}m_{3}^{2}&plus;m_{2}^{2}m_{3}^{2}&space;\right&space;)&plus;K_{2}m_{1}^{2}m_{2}^{2}m_{3}^{2}\right&space;]dV\:&space;\:&space;\:&space;\:&space;\:&space;(3)"><br><br>

となる。この時、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;m_{i}"> は単位磁化ベクトルの成分であり、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;K_{1},&space;K_{2}"> は異方性定数（材料に依存）を表している。<br><br>

## 交換エネルギー<br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{exch}"> は交換エネルギーであり、磁化の向きの空間変動によってのみ決定され、<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{exch}=A\int\left&space;(&space;m_{1,1}^{2}&plus;m_{1,2}^{2}&plus;m_{1,3}^{2}&plus;&space;m_{2,1}^{2}&plus;m_{2,2}^{2}&plus;m_{2,3}^{2}&space;&plus;&space;m_{3,1}^{2}&plus;m_{3,2}^{2}&plus;m_{3,3}^{2}&space;\right&space;)dV\:&space;\:&space;\:&space;\:&space;\:&space;(4)"><br><br>

となる。この時、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;A"> は交換定数である。また、添え字内のコンマは空間微分を表しており、以降
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;m_{i,j}=\frac{\partial&space;m_{i}}{\partial&space;x_{j}}"> とし、
<img src="https://latex.codecogs.com/gif.latex?x_{j}"> をデカルト座標の j 番目の位置ベクトル成分とする。なお、本シミュレーション内では各セルに隣接する上下左右４つのセルに関して計算する。<br><br>

## 静磁エネルギー<br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{ms}"> は静磁エネルギーであり、<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{ms}=-\frac{1}{2}\mu&space;_{0}M_{s}\int&space;\mathbf{H_{d}\cdot&space;m}\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(5)" ><br><br>

となる。この時、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{d}}"> は反磁界であり、<br><br>

<img src="https://latex.codecogs.com/gif.latex?H_{d1,1}&plus;H_{d2,2}&plus;H_{d3,3}=-M_{s}\left&space;(&space;m_{d1,1}&plus;m_{d2,2}&plus;m_{d3,3}&space;\right&space;)\:&space;\:&space;\:&space;\:&space;\:&space;(6)" ><br><br>

と表すことができる。ここで、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{di}}"> は
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{d}}"> の構成要素である（添え字内のコンマは空間微分）。また、<br><br>

<img src="https://latex.codecogs.com/gif.latex?H_{di}=-{\phi&space;}_j&space;\:&space;\:&space;\:&space;\:&space;\:&space;(7)"><br><br>

のように磁気スカラーポテンシャル
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{\phi_{j}}"> を導入することで、(6)式は<br><br>

<img src="https://latex.codecogs.com/gif.latex?\Delta&space;\phi&space;=M_s(m_{1,1}&plus;m_{2,2}&plus;m_{3,3})\:&space;\:&space;\:&space;\:&space;\:&space;(8)" ><br><br>

のように書き直される。この(8)式のポテンシャルの解はフーリエ空間で次のように与えられる。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\phi&space;(k)=-i\frac{M_s\left&space;[&space;m_1(k)k_1&plus;m_2(k)k_2&plus;m_3(k)k_3&space;\right&space;]}{k_{1}^{2}&plus;k_{2}^{2}&plus;k_{3}^{2}}\:&space;\:&space;\:&space;\:&space;\:&space;(9)" ><br><br>

ここで、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;i=\sqrt{-1},&space;\:&space;\:&space;k_i" > はフーリエ空間の座標であり、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\phi(k),&space;\:&space;\:&space;m_i(k)" >  はそれぞれ、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\phi,&space;\:&space;\:&space;m_i" > のフーリエ変換である。つまり実空間での各値は逆フーリエ変換によって求めることができる。これらのことを用いて
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{H_{d}}"> を求めることができる。しかし、式(9)において
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;k_1=k_2=k_3=0"> となる点（フーリエ空間の原点）において無限に発散してまう。この点での寄与は系全体での平均磁化<br><br>

<img src="https://latex.codecogs.com/gif.latex?\mathbf{\bar{M}}=\int\delta&space;M(r)&space;\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(10)">

による反磁界に対応している。つまり、磁化の空間分布は<br><br>

<img src="https://latex.codecogs.com/gif.latex?\mathbf{M(r)}=\mathbf{\bar{M}}+\delta&space;M(\mathbf{r})&space;\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(11)"><br><br>

のように空間に依存する磁化
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\delta&space;M(\mathbf{r})"> と空間に依存しない磁化
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\mathbf{\bar{M}}"> の和で表され、(9)式では空間に依存した磁化の寄与のみを計算している。つまり、
反磁場係数
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;N"> を用いて次の反磁界を追加して考慮する必要がある。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\mathbf{H_d(\bar{M})}=N\mathbf{\bar{M}}\:&space;\:&space;\:&space;\:&space;\:&space;(12)">


## Zeemanエネルギー<br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{extarnal}"> は外部エネルギーであり、外部磁界
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;H_{ex}"> を用いて<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{extarmal}=-\mu_0M_s\int\mathbf{H_{ex}}&space;\cdot&space;\mathbf{m}\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(13)"><br><br>

と表すことができる。<br><br>

## 弾性エネルギー<br>

立方体材料の場合は、局所的な磁化に関する変形は固有ひずみ（外部応力の加わっていないひずみ）によって記述され、<br><br>

<img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;_{11}^{0}=\frac{3}{2}\lambda&space;_{100}\left&space;(&space;m_1^2-\frac{1}{3}&space;\right&space;),\:&space;\:&space;\varepsilon&space;_{12}^{0}=\frac{3}{2}\lambda&space;_{111}m_1m_2"><br>
<img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;_{22}^{0}=\frac{3}{2}\lambda&space;_{100}\left&space;(&space;m_1^2-\frac{1}{3}&space;\right&space;),\:&space;\:&space;\varepsilon&space;_{13}^{0}=\frac{3}{2}\lambda&space;_{111}m_1m_3\:&space;\:&space;\:&space;\:&space;\:&space;(14)"><br>
<img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;_{33}^{0}=\frac{3}{2}\lambda&space;_{100}\left&space;(&space;m_1^2-\frac{1}{3}&space;\right&space;),\:&space;\:&space;\varepsilon&space;_{23}^{0}=\frac{3}{2}\lambda&space;_{111}m_2m_3"><br><br>

となる。ここで、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\lambda_{100}\:&space;,&space;\:&space;\:&space;\lambda_{111}"> は立方晶の磁歪定数である。このとき、磁歪効果から生じる局所的な変形によって弾性ひずみは全ひずみ
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\varepsilon_{ij}"> を用いて、<br><br>

<img src="https://latex.codecogs.com/gif.latex?e_{ij}=\varepsilon&space;_{ij}-\varepsilon&space;_{ij}^0\:&space;\:&space;\:&space;\:&space;\:&space;(15)"><br><br>

と表せられる。またこれに対応する弾性エネルギーは<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{elastic}=\int&space;\frac{1}{2}c_{ijkl}e_{ij}e_{kl}\:&space;dV=\int&space;\frac{1}{2}c_{ijkl}(\varepsilon&space;_{ij}-\varepsilon&space;_{ij}^0)(\varepsilon&space;_{kl}-\varepsilon&space;_{kl}^0)\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(16)"><br><br>

となる。ここで、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;c_{ijkl}"> は弾性係数テンソルである。<br>

ここで弾性係数テンソルに関して説明する。弾性係数テンソルはフックの法則（線形弾性体において）において次のように行列表現を用いて表すことができる。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;\sigma&space;_{11}&space;\\&space;\sigma&space;_{22}&space;\\&space;\sigma&space;_{33}&space;\\&space;\sigma&space;_{23}&space;\\&space;\sigma&space;_{31}&space;\\&space;\sigma&space;_{12}&space;\end{bmatrix}&space;\begin{bmatrix}&space;c_{1111}&space;&c_{1122}&space;&c_{1133}&space;&c_{1123}&space;&c_{1131}&space;&c_{1112}&space;\\&space;(c_{2211})&space;&c_{2222}&space;&c_{2233}&space;&c_{2223}&space;&c_{2231}&space;&c_{2212}&space;\\&space;(c_{3311})&space;&(c_{3322})&space;&c_{3333}&space;&c_{3323}&space;&c_{3331}&space;&c_{3312}&space;\\&space;(c_{2311})&space;&(c_{2322})&space;&(c_{2333})&space;&c_{2323}&space;&c_{2331}&space;&c_{2312}&space;\\&space;(c_{3111})&space;&(c_{3122})&space;&(c_{3133})&space;&(c_{3123})&space;&c_{3131}&space;&c_{3112}&space;\\&space;(c_{1211})&space;&(c_{1222})&space;&(c_{1233})&space;&(c_{1223})&space;&(c_{1231})&space;&c_{1212}&space;\end{bmatrix}&space;\begin{bmatrix}&space;e&space;_{11}&space;\\&space;e&space;_{22}&space;\\&space;e&space;_{33}&space;\\&space;2e&space;_{23}&space;\\&space;2e&space;_{31}&space;\\&space;2e&space;_{12}&space;\end{bmatrix}\:&space;\:&space;\:&space;\:&space;\:&space;(17)"></a>

本来4階のテンソルである弾性係数テンソルは行列の形で表現することはできないが、応力と弾性ひずみの対称性からこのような行列の積で計算することが可能である。また、式(17)内において、弾性係数テンソルは対角行列成分に対して対称的であるので独立な成分は21個となり、Voigt表記を用いると以下のように書き直すことができる。<br><Rb>

<img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;\sigma&space;_{1}&space;\\&space;\sigma&space;_{2}&space;\\&space;\sigma&space;_{3}&space;\\&space;\sigma&space;_{4}&space;\\&space;\sigma&space;_{5}&space;\\&space;\sigma&space;_{6}&space;\end{bmatrix}&space;\begin{bmatrix}&space;c_{11}&space;&c_{12}&space;&c_{13}&space;&c_{14}&space;&c_{15}&space;&c_{16}&space;\\&space;(c_{21})&space;&c_{22}&space;&c_{23}&space;&c_{24}&space;&c_{25}&space;&c_{26}&space;\\&space;(c_{31})&space;&(c_{32})&space;&c_{33}&space;&c_{34}&space;&c_{35}&space;&c_{36}&space;\\&space;(c_{41})&space;&(c_{42})&space;&(c_{43})&space;&c_{44}&space;&c_{45}&space;&c_{46}&space;\\&space;(c_{51})&space;&(c_{52})&space;&(c_{53})&space;&(c_{54})&space;&c_{55}&space;&c_{56}&space;\\&space;(c_{61})&space;&(c_{62})&space;&(c_{63})&space;&(c_{64})&space;&(c_{65})&space;&c_{66}&space;\end{bmatrix}&space;\begin{bmatrix}&space;e&space;_{1}&space;\\&space;e&space;_{2}&space;\\&space;e&space;_{3}&space;\\&space;e&space;_{4}&space;\\&space;e&space;_{5}&space;\\&space;e&space;_{6}&space;\end{bmatrix}\:&space;\:&space;\:&space;\:&space;\:&space;(18)"><br><br>

立方晶においてはその対称性から弾性係数テンソルの独立な成分数は3つとなり、次のように表すことができる。<br><br>

<img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}&space;\sigma&space;_{1}&space;\\&space;\sigma&space;_{2}&space;\\&space;\sigma&space;_{3}&space;\\&space;\sigma&space;_{4}&space;\\&space;\sigma&space;_{5}&space;\\&space;\sigma&space;_{6}&space;\end{bmatrix}&space;\begin{bmatrix}&space;c_{11}&space;&c_{12}&space;&c_{12}&space;&0&space;&0&space;&0&space;\\&space;c_{12}&space;&c_{11}&space;&c_{12}&space;&0&space;&0&space;&0&space;\\&space;c_{12}&space;&c_{12}&space;&c_{33}&space;&0&space;&0&space;&0&space;\\&space;0&space;&0&space;&0&space;&c_{44}&space;&0&space;&0&space;\\&space;0&space;&0&space;&0&space;&0&space;&c_{44}&space;&0&space;\\&space;0&space;&0&space;&0&space;&0&space;&0&space;&c_{44}&space;\end{bmatrix}&space;\begin{bmatrix}&space;e&space;_{1}&space;\\&space;e&space;_{2}&space;\\&space;e&space;_{3}&space;\\&space;e&space;_{4}&space;\\&space;e&space;_{5}&space;\\&space;e&space;_{6}&space;\end{bmatrix}\:&space;\:&space;\:&space;\:&space;\:&space;(19)"><br><br>

式(19)のように３つの成分のみをもつ立方晶の弾性係数テンソルを用いて式(16)の弾性エネルギーを表すと<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{elastic}=\int\left\{\frac{1}{2}c_{11}(&space;e_{11}^{2}&plus;e_{22}^{2}&plus;e_{33}^{2}&space;)&plus;c_{12}(&space;e_{11}e_{22}&plus;e_{22}e_{33}&plus;e_{33}e_{11})&plus;2c_{44}&space;(&space;e_{12}^{2}&plus;e_{23}^{2}&plus;e_{31}^{2})\right\}\:&space;\:&space;\:&space;\:&space;\:&space;(20)"><br><br>

<img src="https://latex.codecogs.com/gif.latex?\dpi{120}&space;\:&space;\:&space;\:&space;\:&space;\:&space;(\:&space;e_{ij}=\varepsilon&space;_{ij}-\varepsilon&space;_{ij}^0)">

となり、以下の３つの寄与に分けることができる。<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{elastic1}=\int\left&space;\{&space;\frac{1}{2}c_{11}(\varepsilon_{11}^{2}&space;&plus;&space;\varepsilon_{22}^{2}&space;&plus;&space;\varepsilon_{33}^{2})&space;&plus;&space;c_{12}(\varepsilon_{11}\varepsilon_{22}&space;&plus;&space;\varepsilon_{22}\varepsilon_{33}&space;&plus;&space;\varepsilon_{33}\varepsilon_{11})&space;&plus;&space;2&space;c_{44}(\varepsilon_{12}^{2}&space;&plus;&space;\varepsilon_{23}^{2}&space;&plus;&space;\varepsilon_{13}^{2})\right&space;\}\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(21)"><br>

<img src="https://latex.codecogs.com/gif.latex?E_{elastic2}=\int\left&space;\{&space;\left&space;[&space;2c_{44}\left&space;(&space;\frac{3}{2}\lambda_{111}&space;\right&space;)^{2}&space;-&space;\left&space;(&space;c_{11}&space;-&space;c{12}&space;\right&space;)\left&space;(&space;\frac{3}{2}\lambda_{100}&space;\right&space;)^{2}&space;\right&space;]\times&space;\left&space;(&space;m_1^2m_2^2&space;&plus;&space;m_2^2m_3^2&space;&plus;&space;m_3^2m_1^2&space;\right&space;)&space;\right&space;\}\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(22)"><br>

<img src="https://latex.codecogs.com/gif.latex?E_{elastic3}=\int\left&space;\{&space;-\frac{3}{2}\lambda_{100}\left&space;(&space;c_{11}&space;-&space;c_{12}&space;\right&space;)\left&space;(&space;\varepsilon_{11}m_1^2&space;&plus;&space;\varepsilon_{22}m_2^2&space;&plus;&space;\varepsilon_{33}m_3^2&space;\right&space;)&space;-6\lambda_{111}c_{44}\left&space;(&space;\varepsilon_{12}m_1m_2&space;&plus;&space;\varepsilon_{23}m_2m_3&space;&plus;&space;\varepsilon_{31}m_3m_1&space;\right&space;)&space;\right&space;\}\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;(23)"><br><br>

Khachatuyanの理論に従うと<br><br>

<img src="https://latex.codecogs.com/gif.latex?\varepsilon_{ij}(\mathbf{r})&space;=&space;\bar{\varepsilon}_{ij}&space;&plus;&space;\eta_{ij}(\mathbf{r})\:&space;\:&space;\:&space;\:&space;\:&space;(24)"><br><br>

のように全ひずみは均一ひずみと不均一ひずみで表され、均一ひずみは<br><br>

<img src="https://latex.codecogs.com/gif.latex?\int\eta_{ij}(\mathbf{r})\:&space;dV&space;=&space;0\:&space;\:&space;\:&space;\:&space;\:&space;(25)"><br><br>

となる（不均一ひずみの体積分が0）ように定義される。つまり均一ひずみは系全体の巨視的な形状変化を表しており、不均一ひずみは巨視的な形状には影響を与えない。平衡不均一ひずみは弾性変位に関してオイラー方程式で与えられる機械的平衡条件<br><br>

<img src="https://latex.codecogs.com/gif.latex?\sigma_{ij,j}&space;=&space;0&space;\Rightarrow&space;\left\{\begin{matrix}&space;\frac{\partial&space;\sigma_{11}&space;}{\partial&space;x_1}&space;&plus;&space;\frac{\partial&space;\sigma_{12}&space;}{\partial&space;x_2}&space;&plus;&space;\frac{\partial&space;\sigma_{31}&space;}{\partial&space;x_3}&space;=&space;0&space;\\&space;\frac{\partial&space;\sigma_{21}&space;}{\partial&space;x_1}&space;&plus;&space;\frac{\partial&space;\sigma_{22}&space;}{\partial&space;x_2}&space;&plus;&space;\frac{\partial&space;\sigma_{23}&space;}{\partial&space;x_3}&space;=&space;0\\&space;\frac{\partial&space;\sigma_{31}&space;}{\partial&space;x_1}&space;&plus;&space;\frac{\partial&space;\sigma_{32}&space;}{\partial&space;x_2}&space;&plus;&space;\frac{\partial&space;\sigma_{33}&space;}{\partial&space;x_3}&space;=&space;0&space;\end{matrix}\right.\:&space;\:&space;\:&space;\:&space;\:&space;(26)" ><br><br>

を満たし、応力成分
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\sigma_{ij}"> は(17)の行列式で与えられる。<br>
均一な弾性係数で近似を行う場合、平衡不均一ひずみはフーリエ空間で(24)式を解くことによって求めることができる。まず、変位
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;u_{i}(\mathbf{r})"> のセットを導入すると、<br><br>

<img src="https://latex.codecogs.com/gif.latex?\eta_{ij}&space;=&space;\frac{1}{2}\left&space;(&space;u_{i,j}&space;&plus;&space;u_{j,i}&space;\right&space;)\:&space;\:&space;\:&space;\:&space;\:&space;(27)" ><br><br>










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

エネルギー項は順に結晶磁気異方性、交換、静磁、外部、弾性エネルギーであり、下記に求め方を順に示す。また、このプログラム内では簡単のため、 
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


## 外部エネルギー<br>

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;E_{extarnal}"> は外部エネルギーであり、外部磁界
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;H_{ex}"> を用いて<br><br>

<img src="https://latex.codecogs.com/gif.latex?E_{extarmal}=-\mu_0M_s\int\mathbf{H_{ex}}&space;\cdot&space;\mathbf{m}\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;()"><br><br>

と表すことができる。<br><br>

## 弾性エネルギー<br>

立方体材料の場合は、局所的な磁化に関する変形は固有ひずみ（応力を常時ないひずみ）によって記述され、<br><br>

<img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;_{11}^{0}=\frac{3}{2}\lambda&space;_{100}\left&space;(&space;m_1^2-\frac{1}{3}&space;\right&space;),\:&space;\:&space;\varepsilon&space;_{12}^{0}=\frac{3}{2}\lambda&space;_{111}m_1m_2"><br>
<img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;_{22}^{0}=\frac{3}{2}\lambda&space;_{100}\left&space;(&space;m_1^2-\frac{1}{3}&space;\right&space;),\:&space;\:&space;\varepsilon&space;_{13}^{0}=\frac{3}{2}\lambda&space;_{111}m_1m_3\:&space;\:&space;\:&space;\:&space;\:&space;()"><br>
<img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;_{33}^{0}=\frac{3}{2}\lambda&space;_{100}\left&space;(&space;m_1^2-\frac{1}{3}&space;\right&space;),\:&space;\:&space;\varepsilon&space;_{23}^{0}=\frac{3}{2}\lambda&space;_{111}m_2m_3"><br><br>

となる。ここで、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;\lambda_{100}\:&space;,&space;\:&space;\:&space;\lambda_{111}"> は立方晶の磁歪定数である。このとき、磁歪効果から生じる局所的な変形によって弾性ひずみは全ひずみ
 を用いて、<br><br>

<img src="https://latex.codecogs.com/gif.latex?e_{ij}=\varepsilon&space;_{ij}-\varepsilon&space;_{ij}^0\:&space;\:&space;\:&space;\:&space;\:&space;()"><br><br>

と表せられる。またこれに対応する弾性エネルギーは

<img src="https://latex.codecogs.com/gif.latex?E_{elastic}=\int&space;\frac{1}{2}c_{ijkl}e_{ij}e{kl}\:&space;dV=\int&space;\frac{1}{2}c_{ijkl}(\varepsilon&space;_{ij}-\varepsilon&space;_{ij}^0)(\varepsilon&space;_{kl}-\varepsilon&space;_{kl}^0)\:&space;dV\:&space;\:&space;\:&space;\:&space;\:&space;()
">

となる。ここで、
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{80}&space;c_{ijkl}"> は弾性係数テンソルであり、次のように行列表現を用いて表すことができる。


线性回归

模型：\mathbit{y}=\mathbit{Xw}+\mathbit{\varepsilon}, 即X\in R_{n\times k}，\mathbit{w}\in R_{k\times1}，\mathbit{\varepsilon}\in R_{n\times1}\ （有n个样本，k个变量）

五条基本假设：
1，存在线性关系y=Xw+\varepsilon	
2，X是确定值，且满秩可逆，即不存在共线性
3，误差项均值为0，E\left(\varepsilon\right)=0
4，误差项同方差，Var(\varepsilon_1)\ =\ Var(\varepsilon_2)\ =\ \cdots\ =\ Var(\varepsilon_n)\ =\ \sigma^2
5, 误差项不相关，即Cov(\varepsilon_i)\ \neq\ Cov(\varepsilon_j)\ 
3,4,5条可以改写为\varepsilon服从N(0,\sigma^2I_n)，即\varepsilon/\sigma服从N(0,I_n)

计算估计值

估计值 \hat{\omega}, \hat{\varepsilon}。 通过最小二乘法求解\hat{\omega}，通过y=X\hat{\omega}+\hat{\varepsilon}求解\hat{\varepsilon}
\hat{\mathbit{\omega}}\ =\ {{(\mathbit{X}}^\mathbit{T}\mathbit{X})}^{-\mathbf{1}}\mathbit{X}^\mathbit{T}\mathbit{y}\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ （1）
	E\left(\hat{\omega}\right)=w, \hat{\omega}是w的无偏估计
	Var(\hat{\omega})\ =\ \sigma^2{{(X}^TX)}^{-1}
{\hat{\mathbit{\omega}}}_\mathbit{j}~\mathbit{N}(\mathbit{w}_\mathbit{j},\ \sigma^2{{{(X}^TX)}^{-1}}_{jj})\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (2)
\hat{\varepsilon}=\ y\ -\ X\hat{\omega}\ =\ y\ -\ X{{(X}^TX)}^{-1}X^Ty\ =\ (I\ -\ X{{(X}^TX)}^{-1}X^T)\ y\ =\ My\ 
这里M=\ I\ -\ X{{(X}^TX)}^{-1}X^T 
易证\ MX=\ (I\ -\ X{{(X}^TX)}^{-1}X^T)X\ =\ 0
\hat{\mathbit{\varepsilon}}=\mathbit{My}=\mathbit{M}\left(\mathbit{Xw}+\mathbit{\varepsilon}\right)=\mathbit{M\varepsilon}\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 3
易证 \mathbit{M}是幂等矩阵，对称矩阵
M^T=\ {(I\ -\ X{{(X}^TX)}^{-1}X^T)}^T\ =\ I\ -\ X{{(X}^TX)}^{-1}X^T\ =\ M
MM=\left(I-X{{(X}^TX)}^{-1}X^T\right)\left(I-X{{(X}^TX)}^{-1}X^T\right)=I-X{{(X}^TX)}^{-1}X^T=\ M\ 
{\hat{\mathbit{\varepsilon}}}^\mathbit{T}\hat{\mathbit{\varepsilon}}\ =\ \mathbit{tr}({\hat{\mathbit{\varepsilon}}}^\mathbit{T}\hat{\mathbit{\varepsilon}}), 由于{\hat{\mathbit{\varepsilon}}}^\mathbit{T}\hat{\mathbit{\varepsilon}}是一个标量。另由（3）式不难得出 {\hat{\mathbit{\varepsilon}}}^\mathbit{T}\hat{\mathbit{\varepsilon}}={\ \mathbit{\varepsilon}}^\mathbit{T}\mathbit{M}^\mathbit{T}\mathbit{M\varepsilon}=\ {\ \mathbit{\varepsilon}}^\mathbit{T}\mathbit{M\varepsilon}\ 
离差平方和
tr(M)\ =\ tr(I_n-\ X{{(X}^TX)}^{-1}X^T)=tr(I_n)-tr(X{{(X}^TX)}^{-1}X^T)\ =\ n\ -\ tr(X^TX{{(X}^TX)}^{-1})=\ n-\ tr(I_k)\ =\ n-k
E(RSS)\ =\ E({\hat{\varepsilon}}^T\hat{\varepsilon})={\ E(\varepsilon}^TM^TM\varepsilon)=\ {\ E(\varepsilon}^TM\varepsilon)
E(RSS)\ =\ E({\hat{\varepsilon}}^T\hat{\varepsilon})=E(tr({\hat{\varepsilon}}^T\hat{\varepsilon}))=\ {\ E(tr(\varepsilon}^TM\varepsilon))=\ E(tr(M{\ \varepsilon}^T\varepsilon\ )=E(tr(M\sigma^2I)=\sigma^2E(tr(M))\ =\ (n-k)\sigma^2
注：矩阵的迹（对角线各元素的和），trace()或者tr()，重要两条性质
tr\left(AB\right)=tr\left(BA\right), tr\left(A+B\right)=tr\left(A\right)+tr\left(B\right)
即E({\hat{\sigma}}^2)=E\left({\hat{\varepsilon}}^T\hat{\varepsilon}/(n-k)\right)=\sigma^2, {\hat{\sigma}}^2={\hat{\varepsilon}}^T\hat{\varepsilon}/(n-k)是\sigma^2的无偏估计
分布
\frac{\varepsilon}{\sigma}~N(0,I_n) => \frac{\varepsilon^T}{\sigma}M\frac{\varepsilon}{\sigma}\ ~\chi^2(n-k)
证明：令z=\frac{\varepsilon}{\sigma}，即\frac{\varepsilon^T}{\sigma}M\frac{\varepsilon}{\sigma}\ =\ z^TMz；由M是对称矩阵，即存在正交矩阵C，满足正交变换M=\ C\Lambda C^T,
\Lambda是对角矩阵，对角线上元素是特征值，由于M是幂等矩阵，即M^n\ =\ C\mathrm{\Lambda}^nC^T=\ C\Lambda C^T, 所有对角线上元素只能为0或者1。\frac{\varepsilon^T}{\sigma}M\frac{\varepsilon}{\sigma}\ =\ z^TMz=\ z^TC\Lambda C^Tz, 由正态分布的性质，x\ ~\ N(\mu,\sum)， 则Ax+b\ ~\ N(A\mu+b,A\sum A^T)
z=\frac{\varepsilon}{\sigma}~N(0,I_n)=> z=\frac{\varepsilon}{\sigma}~N(0,I_n) => C^Tz\ ~\ N(0,C^TI_nC)\ =N(0,I_n), \frac{\varepsilon^T}{\sigma}M\frac{\varepsilon}{\sigma}是n-k个标准正态分布相加，服从\chi^2(n-k)
\frac{\varepsilon^T}{\sigma}M\frac{\varepsilon}{\sigma},\ \frac{{\hat{\varepsilon}}^T\hat{\varepsilon}}{\sigma^2},\ \frac{RSS}{\sigma^2}\ ,\ \frac{(n-k){\hat{\sigma}}^2}{\sigma^2}~\chi^2(n-k)\ \ \ \ \ \ \ \ (4)
T检验
由（2）可知\frac{{\hat{\mathbit{\omega}}}_\mathbit{j}-\mathbit{w}_\mathbit{j}}{\sqrt{\sigma^2{{{(X}^TX)}^{-1}}_{jj}}}~\mathbit{N}\left(\mathbf{0},\mathbf{1}\right)，\sfrac{\frac{{\hat{\omega}}_j-w_j}{\sqrt{\sigma^2{{{(X}^TX)}^{-1}}_{jj}}}}{\sqrt{\sfrac{\frac{\left(n-k\right){\hat{\sigma}}^2}{\sigma^2}}{\left(n-k\right)}}}~t\left(n-k\right)
=> \frac{{\hat{\mathbit{\omega}}}_\mathbit{j}-\mathbit{w}_\mathbit{j}}{\sqrt{{{{{\hat{\mathbit{\sigma}}}^\mathbf{2}(\mathbit{X}}^\mathbit{T}\mathbit{X})}^{-\mathbf{1}}}_{\mathbit{jj}}}}\ ~\mathbit{t}\left(\mathbit{n}-\mathbit{k}\right)， {{{{\hat{\sigma}}^2(X}^TX)}^{-1}}_{jj}是w_j方差{Var(w}_j) 的估计值\widehat{{Var(w}_j)}，分母称为w_j的标准误standard error(s.e.).
因此t_j=\frac{{\hat{\mathbit{\omega}}}_\mathbit{j}}{\mathbit{s}.\mathbit{e}.(w_j)}用来检验w_j是否为0
F检验
上面t检验是检验单个w_j，若检验多个w_j，如y=w_0+w_1x_1+w_2x_2+w_3x_3，检验模型整体显著，H_0:\ w_1=\ w_2=w_3\ =0，可用Rw-r=0可写\left[\begin{matrix}0\\0\\0\\\end{matrix}\begin{matrix}1\\0\\0\\\end{matrix}\begin{matrix}0\\1\\0\\\end{matrix}\begin{matrix}0\\0\\1\\\end{matrix}\right]\left[\begin{matrix}w_0\\w_1\\w_2\\w_3\\\end{matrix}\right]-\left[\begin{matrix}0\\0\\0\\\end{matrix}\right]=\left[\begin{matrix}0\\0\\0\\\end{matrix}\right]，这里R是p\times k矩阵，r是p\times1向量。
由前面标绿色的证明易知，对于p\times1向量g~N(\mu,\mathrm{\Sigma})，
g-μ)TΣ-1(g-μ ~ χ2(p)
R\hat{w}-r\ ~\ N(Rw-r,\ R\sigma^2{{(X}^TX)}^{-1}R^T)
由\ Rw-r=0，Rw-rT[Rσ2(XTX)-1RT]-1Rw-r~χ2(p)
即Rw-rT[R(XTX)-1RT]-1Rw-rσ2~χ2p
由（4）式\frac{(n-k){\hat{\sigma}}^2}{\sigma^2}~\chi^2(n-k)，及F分布定义（若随机变量X~\chi^2(a), Y~\chi^2(b)，则 \frac{\sfrac{X}{a}}{\sfrac{Y}{b}}\ ~\ F_{a,b}）得 Rw-rT[R(XTX)-1RT]-1Rw-rσ2pn-kσ2σ2n-k=Rw-rT[R(XTX)-1RT]-1Rw-rσ2p
由{\hat{\sigma}}^2=RSS/(n-k)\ =\ {MS}_E，Rw-rT[R(XTX)-1RT]-1Rw-rp= SSRp=MSR, \frac{{MS}_R}{{MS}_E} ~ F_{p,\ n-k}, p其实就是字变量个数k。
SSR


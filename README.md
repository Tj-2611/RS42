# RS42
The repository contains the codes in the manuscript "4/2 Rough and Smooth" (see [https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4909939])
We sincerely thank the authors of the following works and their codes, which have greatly helped this research.
> 1. Aït-Sahalia, Yacine, and Robert Kimmel. "Maximum likelihood estimation of stochastic volatility models." Journal of Financial Economics 83, no. 2 (2007): 413-452.
See [https://www.princeton.edu/~yacine/research]
> 2. Baschetti, Fabio, Giacomo Bormetti, Silvia Romagnoli, and Pietro Rossi. "The SINC way: A fast and accurate approach to Fourier pricing." Quantitative Finance 22, no. 3 (2022): 427-446.
See [https://github.com/fabioBaschetti/SINC-method]
> 3. Markus Buehren. "Differential Evolution" (2019). MATLAB Central File Exchange
See [https://www.mathworks.com/matlabcentral/fileexchange/18593-differential-evolution]
> 4. Gatheral, Jim, and Radoš Radoičić. "A generalization of the rational rough Heston approximation." Quantitative Finance 24, no. 2 (2024): 329-335.
See [https://github.com/jgatheral/rationalroughheston]
> 5. Rømer, Sigurd Emil. "Empirical analysis of rough and classical stochastic volatility models to the SPX and VIX markets." Quantitative Finance 22, no. 10 (2022): 1805-1838.
See [https://github.com/sigurdroemer/rough_volatility]

## To run the codes:
1. You may first download the volatility surfaces data from OptionMetrics for SPX and VIX options and concatenate them to obtain two files:
"surf_filt_liq.csv" and "vix_option_future.csv". You may move them into the directory '/data/' 
See https://wrds-www.wharton.upenn.edu/pages/about/data-vendors/optionmetrics/ 
For SPX options, each quote is matched with the closing prices of S\&P500 index, which is adjusted for dividends using the dividend rates
in OptionMetrics. For VIX options, we consider VIX futures to be the underlying level whose prices are
inferred from highly liquid options using ATM put–call parity. The risk-free interest rate for each option
maturity is calculated by interpolating the zero-coupon interest rate curve from OptionMetrics. We also
download the sensitivities, namely Vega, of the filtered SPX and VIX options from OptionMetrics.
See our paper for more details. We have added two sample files in '/data/' just for reference (there is no real data).

2. Download "Differential Evolution" by Markus Buehren. Unzip it and rename it as 'DE' in the directory of our project. 
See [https://ww2.mathworks.cn/matlabcentral/fileexchange/18593-differential-evolution]
 
## What can these codes do:
1. in_sample: two-step iterative estimation for 4/2RS model, 4/2SV model, H2F model. You will obtain the estimates results in a directory '/estimates/';
2. out_of_sample: out-of-sample test based on the structural parameters estimated in 1;
3. analysis: some related analyses such as risk-neutral log-return moment estimation, variance risk premium estimation, model-implied EVs based on simulations; 
4. source: key intermediate codes in this project.



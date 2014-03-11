********************************************************************************
** Modeling Cyclical Patterns in Daily College Drinking Data with Many Zeroes **
**                                                                            **
** Authors: Huh, D., Kaysen, D., & Atkins, D.C.                               **
** Software: SAS 9.4                                                          **
*******************************************************************************;

options formdlim=' ' ls=90 ps=50 nofmterr;


********** Import Data **************************;
filename dashweb url 'https://raw.github.com/davidhuh/cyclicalmodels/master/dash.csv';
proc import datafile=dashweb out=WORK.dash dbms=CSV;
     getnames=YES;
     datarow=2;
run;


********** Prepare Data for Analysis ************;
data dash; set dash;
  ** Outcome variable
  y = drinks;
  
  ** Create cyclical covariates;
  pi = constant('PI');
  
  cycPha = sin(2*pi*(dayofwk/7));
  cycAmp = pha(2*pi*(dayofwk/7));
run;

** Sort observations by participant ID;
proc sort data=dash;
  by id;
run;

** Send SAS output to PDF file;
ods pdf file='H:\CSHRB\Cyclical Models\HurdleNB_SAS.pdf';


********** Hurdle Negative Binomial Mixed Model -- No Time Predictors ******;

** Model without random effects for starting values ****;
proc nlmixed data=dash;
  title 'Hurdle Negative Binomial (NB2) Fixed Effect Model';
  parms a0=1 a1=0 a2=0 a3=0
        b0=1 b1=0 b2=0 b3=0;

  ** Logit model;
  eta_zip= a0 + a1*dmqsoc + a2*cycPha + a3*cycAmp;
  p0_zipe= 1/(1+exp(-1*eta_zip));

  ** Truncated count model;
  eta_nb = b0 + b1*dmqsoc + b2*cycPha + b3*cycAmp;
  mean=exp(eta_nb);
  p=1/(1+(1/v)*mean);

  p0=log(p0_zipe);
  p_else=log(1-p0_zipe)+y*log(1-p)-log(p**(-1*(v))-1)+lgamma(y+(v))-
  lgamma(v)-log(fact(y));

  if y=0 then loglike=(p0);
  else loglike=(p_else);

  model y ~ general(loglike);
  logtheta = log(v);

  ods output ParameterEstimates=peHUNB_init1;

quit;title;

** FINAL model with Random Effects ****;
proc nlmixed data=dash;

  title 'Hurdle Negative Binomial (NB2) Mixed Effect Model';
  bounds cv11 cv22 cv33 cv44 > 0;
  parms /data=peHUNB_init1;

  ** Logit model -- Random intercepts;
  eta_zip= a0 + a1*dmqsoc + a2*cycPha + a3*cycAmp + u1;
  p0_zipe= 1/(1+exp(-1*eta_zip));

  ** Truncated count model -- Random intercepts and slopes;
  eta_nb = b0 + b1*dmqsoc + b2*cycPha + b3*cycAmp + u2 + u3*cycPha + u4*cycAmp;
  mean=exp(eta_nb);
  p=1/(1+(1/v)*mean);

  p0=log(p0_zipe);
  p_else=log(1-p0_zipe)+y*log(1-p)-log(p**(-1*(v))-1)+lgamma(y+(v))-
  lgamma(v)-log(fact(y));

  if y=0 then loglike=(p0);

  else loglike=(p_else);

  model y ~ general(loglike);
  logtheta = log(v);

  ** Diagonal covariance matrix of random effects;
  random u1 u2 u3 u4 ~ normal([0,0,0,0],[cv11,0,cv22,0,0,cv33,0,0,0,cv44]) subject=ID;

quit;title;


********** Hurdle Negative Binomial Mixed Effect Model - w/ Time Predictor ********;

** Model without random effects for starting values ****;
proc nlmixed data=dash;
  title 'Hurdle Negative Binomial (NB2) Fixed Effect Model';
  parms a0=1 a1=0 a2=0 a3=0 a4=0 a5=0
        b0=1 b1=0 b2=0 b3=0 b4=0 b5=0;

  ** Logit model;
  eta_zip= a0 + a1*dmqsoc + a2*cycPha + a3*cycAmp +
           a4*cycPha*dmqsoc + a5*cycAmp*dmqsoc;
  p0_zipe=1/(1+exp(-1*eta_zip));

  ** Truncated count model;
  eta_nb = b0 + b1*dmqsoc + b2*cycPha + b3*cycAmp +
           b4*cycPha*dmqsoc + b5*cycAmp*dmqsoc;
  mean=exp(eta_nb);
  p=1/(1+(1/v)*mean);

  p0=log(p0_zipe);
  p_else=log(1-p0_zipe)+y*log(1-p)-log(p**(-1*(v))-1)+lgamma(y+(v))-
  lgamma(v)-log(fact(y));

  if y=0 then loglike=(p0);
  else loglike=(p_else);

  model y ~ general(loglike);
  logtheta = log(v);

  ods output ParameterEstimates=peHUNB_init2;

quit;title;

** FINAL model with Random Effects ****;
proc nlmixed data=dash;

  title 'Hurdle Negative Binomial (NB2) Mixed Effect Model';
  bounds cv11 cv22 cv33 cv44 > 0;
  parms /data=peHUNB_init2;

  ** Logit model -- Random intercepts;
  eta_zip= a0 + a1*dmqsoc + a2*cycPha + a3*cycAmp +
           a4*cycPha*dmqsoc + a5*cycAmp*dmqsoc + u1;
  p0_zipe=1/(1+exp(-1*eta_zip));

  ** Truncated count model -- Random intercepts and slopes;
  eta_nb = b0 + b1*dmqsoc + b2*cycPha + b3*cycAmp +
           b4*cycPha*dmqsoc + b5*cycAmp*dmqsoc + u2 + u3*cycPha + u4*cycAmp;
  mean=exp(eta_nb);
  p=1/(1+(1/v)*mean);

  p0=log(p0_zipe);
  p_else=log(1-p0_zipe)+y*log(1-p)-log(p**(-1*(v))-1)+lgamma(y+(v))-
  lgamma(v)-log(fact(y));

  if y=0 then loglike=(p0);

  else loglike=(p_else);

  model y ~ general(loglike);
  logtheta = log(v);

  ** Diagonal covariance matrix of random effects;
  random u1 u2 u3 u4 ~ normal([0,0,0,0],[cv11,0,cv22,0,0,cv33,0,0,0,cv44]) subject=ID;

quit;title;

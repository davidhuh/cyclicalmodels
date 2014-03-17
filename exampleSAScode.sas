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
  ** Outcome variable;
  y = drinks;
  
  ** Create cyclical covariates;
  pi = constant('PI');
  
  cycPha = sin(2*pi*(dayofwk/7));
  cycAmp = cos(2*pi*(dayofwk/7));

  ** Create saturated set of dummy variables;
  Tue = (dayofwk=1);
  Wed = (dayofwk=2);
  Thu = (dayofwk=3);
  Fri = (dayofwk=4);
  Sat = (dayofwk=5);
  Sun = (dayofwk=6);
run;

** Sort observations by participant ID;
proc sort data=dash;
  by id;
run;

** Send SAS output to PDF file;
ods pdf file='HurdleNB_SAS.pdf';


********** Hurdle Negative Binomial Mixed Model -- No Time Predictors ******;

** Model without random effects for starting values ****;
proc nlmixed data=dash;
  title 'Hurdle Negative Binomial (NB2) Fixed Effect Model';
  parms a0=1 a1=0 a2=0 a3=0 a4=0 a5=0 a6=0 
        b0=1 b1=0 b2=0 b3=0;

  ** Logit model;
  eta_zip= a0 + a1*Tue + a2*Wed + a3*Thu + a4*Fri + a5*Sat + a6*Sun;
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
  eta_zip= a0 + a1*Tue + a2*Wed + a3*Thu + a4*Fri + a5*Sat + a6*Sun + u1;
  p0_zipe= 1/(1+exp(-1*eta_zip));

  ** Truncated count model -- Random intercepts and slopes;
  eta_nb = b0 + b1*dmqsoc + b2*cycAmp + b3*cycPha + u2 + u3*cycAmp + u4*cycPha;
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

  ** Logistic submodel w/ Saturated set of dummy variables;
  **  -- Random intercepts;
  eta_zip= a0 + a1*dmqsoc + a2*Tue + a3*Wed + a4*Thu + a5*Fri + a6*Sat + a7*Sun +
           a8*Tue*dmqsoc + a9*Wed*dmqsoc + a10*Thu*dmqsoc +
           a11*Fri*dmqsoc + a12*Sat*dmqsoc + a13*Sun*dmqsoc;
  p0_zipe=1/(1+exp(-1*eta_zip));

  ** Truncated Negative Binomial (NB2) submodel w/ Cyclical variables;
  **  -- Random intercepts and slopes;
  eta_nb = b0 + b1*dmqsoc + b2*cycAmp + b3*cycPha +
           b4*cycAmp*dmqsoc + b5*cycPha*dmqsoc;
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

  ** Logistic submodel w/ Saturated set of dummy variables;
  eta_zip= a0 + a1*dmqsoc + a2*Tue + a3*Wed + a4*Thu + a5*Fri + a6*Sat + a7*Sun +
           a8*Tue*dmqsoc + a9*Wed*dmqsoc + a10*Thu*dmqsoc +
           a11*Fri*dmqsoc + a12*Sat*dmqsoc + a13*Sun*dmqsoc + u1;
  p0_zipe=1/(1+exp(-1*eta_zip));

  ** Truncated Negative Binomial (NB2) submodel w/ Cyclical variables;
  eta_nb = b0 + b1*dmqsoc + b2*cycAmp + b3*cycPha +
           b4*cycAmp*dmqsoc + b5*cycPha*dmqsoc + u2 + u3*cycAmp + u4*cycPha;
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

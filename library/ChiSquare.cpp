// 自由度nのカイ2乗の上側確率
double qchi(double x2, int n)
{
    double w, pw, x, qc;
    int i, i1;

    if(n < 1) {
        fprintf(stderr,"Error : 自由度 < 1 in qchi()!\n");
        return 0.;
    }
    if(x2 <= 0.)    return 1.;
    if(x2 > 400.)   return 0.;
    if(n > 10) {
        w = 2. / (9. * (double)n);
        pw = pow(x2 / (double)n, 1. / 3.);
        return 1. - qnorm((pw - 1. + w) / sqrt(w));
    }
    w = exp(-x2 / 2.);
    if(n == 2)  return w;
    x = sqrt(x2);
    if(n == 1)  return 2. * (1. - qnorm(x));
    if((n % 2) == 0) {
        i1 = 2;
        qc = w;
    } else {
        i1 = 1;
        qc = 2. * (1. - qnorm(x));
        w *= (0.797884560750774 / x);
    }
    for(i = i1; i <= n - 2; i += 2) {
        w *= (x2 / (double)i);
        qc += w;
    }
    return qc;
}

// 自由度nのカイ2乗のパーセント点
double pchi(double qc, int n)
{
    double c1, c2, gam, x, w, wx;
    int i, j;

    if(qc <= 0. || qc >= 1. || n < 1) {
        fprintf(stderr,"Error : Illigal parameter in pchi()!\n");
        return 0.;
    }
    if(n == 1) {
        w = pnorm(qc / 2.);
        return w * w;
    }
    if(n == 2)  return -2. * log(qc);

    x = -pnorm(qc);
    if(n > 10) {
        w = x * x;
        wx = sqrt(2. * (double)n);
        c1 = (w - 7.) * x / 9. / wx;
        c2 = ((3. * w + 7.) * w - 16.) * 2. / 405. / (double)n;
        wx = (double)n + wx * x + 2. * (w - 1.) / 3. + c1 - c2;
        if(wx < 0.) return 0.;
        return wx;
    }

    w = 2. / 9. / (double)n;
    w = 1. - w + x * sqrt(w);
    wx = (double)n * w * w * w;
    if((n % 2) == 0)    gam = 1.;
    else                gam = 1.772453850905516;
    j = (n + 1) / 2 - 1;
    w = (double)n / 2.;
    for(i = 1; i <= j; i++) gam *= (w - (double)i);
    x = wx / 2.;
    c1 = pow(x, w - 1.);
    c2 = exp(-x) / 2.;
    return wx + (qchi(wx, n) - qc) * gam / c1 / c2;
}



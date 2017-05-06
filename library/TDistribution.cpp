// 正規分布のパーセント点
// 戸田の近似式
double pnorm(double qn)
{
    static double b[11] = {  1.570796288,     0.03706987906,  -0.8364353589e-3,
        -0.2250947176e-3, 0.6841218299e-5, 0.5824238515e-5,
        -0.104527497e-5,  0.8360937017e-7,-0.3231081277e-8,
        0.3657763036e-10,0.6936233982e-12};
    double w1, w3;
    int i;

    if(qn < 0. || 1. < qn) {
        fprintf(stderr, "Error : qn <= 0 or qn >= 1  in pnorm()!\n");
        return 0.;
    }
    if(qn == 0.5)   return 0.;

    w1 = qn;
    if(qn > 0.5)    w1 = 1. - w1;
    w3 = -log(4. * w1 * (1. - w1));
    w1 = b[0];
    for(i = 1; i < 11; i++) w1 += (b[i] * pow(w3, (double)i));
    if(qn > 0.5)    return sqrt(w1 * w3);
    return -sqrt(w1 * w3);
}
// 自由度nのt分布の両側確率
void qtsub(double *q, int n, double w2, double w3, double t2)
{
    int i, j;
    double w;

    j = (n - 2) / 2;
    for(i = 1; i <= j; i++) {
        w = 2. * (double)i - w3;
        *q -= (w2 *= (t2 * w / (1. + w)));
    }
}
double qt(double t, int n)
{
    double t1, t2, w1, w2, wq;

    if(n < 1)
    {
        fprintf(stderr,"Error : n < 1  in qt()!\n");
        return 0.;
    }

    w1 = 0.636619772284456;
    if(t < 0.)  t = - t;
    t1 = t / sqrt((double)n);
    t2 = 1. / (1. + t1 * t1);
    if((n % 2) != 0) {
        wq = 1. - w1 * atan(t1);
        if(n != 1) {
            wq -= (w2 = w1 * t1 * t2);
            if (n != 3) qtsub(&wq, n, w2, 0., t2);
        }
        if(wq > 0.) return wq;
        return 0.;
    }
    wq = 1. - (w2 = t1 * sqrt(t2));
    if(n != 2)      qtsub(&wq, n, w2, 1., t2);
    if(wq > 0.) return wq;
    return 0.;
}


// 自由度nのt分布のパーセント点
double ptsub(double q, int n)
{
    double eps, qe, s, w;

    if(n == 1 && 0.001 <= q && q < 0.01)    eps = 1.e-4;
    else if (n == 2 && q < 0.0001)          eps = 1.e-4;
    else if (n == 1 && q < 0.001)           eps = 1.e-2;
    else                                    eps = 1.e-5;
    s = 10000.;
    w = 0.;
    for(;;) {
        w += s;
        if(s <= eps)    return w;
        if((qe = qt(w,n) - q) == 0.)    return w;
        if(qe < 0.) {
            w -= s;
            s /= 10.;
        }
    }
} 
double pt(double q, int n)
{
    double f, f1, f2, f3, f4, f5, u, u2, w, w0, w1, w2, w3, w4;

    if(q < 1.e-5 || q > 1. || n < 1) {
        fprintf(stderr,"Error : Illigal parameter  in pt()!\n");
        return 0.;
    }

    if(n <= 5)  return ptsub(q, n);

    if(q <= 5.e-3 && n <= 13)   return ptsub(q, n);

    f1 = 4. * (f = (double)n);
    f5 = (f4 = (f3 = (f2 = f * f) * f) * f) * f;
    f2 *= 96.;
    f3 *= 384.;
    f4 *= 92160.;
    f5 *= 368640.;
    u = pnorm(1. - q / 2.);

    w0 = (u2 = u * u) * u;
    w1 = w0 * u2;
    w2 = w1 * u2;
    w3 = w2 * u2;
    w4 = w3 * u2;
    w = ((w0 + u) / f1);
    w += ((5. * w1 + 16. * w0 + 3. * u) / f2);
    w += ((3. * w2 + 19. * w1 + 17. * w0 - 15. * u) / f3);
    w += ((79. * w3 + 776. * w2 + 1482. * w1 - 1920. * w0 - 945. * u) / f4);
    w += ((27. * w4 + 339. * w3 + 930. * w2 - 1782. * w1 - 765. * w0
                + 17955. * u) / f5);
    return u + w;
}



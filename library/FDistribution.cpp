// 自由度n1, n2のF分布の上側確率
double qf(double f, int n1, int n2)
{
    double d, p, w, y, z;
    int j, k, m, mn, n;

    if(f < 0. || n1 < 1 || n2 < 1) {
        fprintf(stderr,"Error : Illegal parameter  in qf()!\n");
        return 0.;
    }

    w = f * (double)n1 / (double)n2;
    z = 1. / (1. + w);
    m = 2 * (n1 / 2) - n1 + 2;
    n = 2 * (n2 / 2) - n2 + 2;
    if((mn = 4 - 2 * (n1 % 2) - (n2 % 2)) == 1) {
        y = 0.318309886142228;
        p = sqrt(w);
        d = y * z / p;
        p = 2. * y * atan(p);
    } else if(mn == 2) {
        p = sqrt(z * w);
        d = p * z / w / 2.;
    } else if(mn == 3) {
        p = sqrt(z);
        d = z * p / 2.;
        p = 1. - p;
    } else {
        d = z * z;
        p = z * w;
    }
    if(n1 <= 2 && n2 <= 2)  return 1. - p;

    y = 2. * w / z;
    if((k = n + 2) <= n2) {
        for( ; k <= n2; k += 2) {
            d *= ((1. + (double)m / (double)(k - 2)) * z);
            if(m != 1)  p = (p + w) * z;
            else        p += (d * y / (double)(k - 1));
        }
    }
    k = m + 2;
    if(k > n1)  return 1. - p;

    y = z * w;
    z = 2. / z;
    j = n2 - 2;
    for( ; k <= n1; k += 2) {
        d *= (y * (double)(k + j) / (double)(k - 2));
        p -= (z * d / (double)(k + j));
    }
    return 1. - p;
}

double pfsub(double x, double y, double z) {
    return (sqrt(z) - y) / x / 2.;
}

// 自由度n1, n2のF分布のパーセント点
double pf(double q, int n1, int n2)
{
    double a, b, c, d, eps, fw, qe, qn, s, u, u2, w1, w2, w3, w4;

    if(q < 0. || q > 1. || n1 < 1 || n2 < 1) {
        fprintf(stderr,"Error : Illegal parameter  in pf()!\n");
        return 0.;
    }

    if(n1 <= 240 || n2 <= 240) {
        eps = 1.e-5;
        if(n2 == 1) eps = 1.e-4;
        s = 1000.;
        fw = 0.;
        for(;;) {
            fw += s;
            if(s <= eps)    return fw;
            if((qe = qf(fw,n1,n2) - q) == 0.)   return fw;
            if(qe < 0.) {
                fw -= s;
                s /= 10.;
            }
        }
    }

    eps = 1.e-6;
    qn = q;
    if(q < 0.5) qn = 1. - q;
    u = pnorm(qn);
    w1 = 2. / (double)n1 / 9.;
    w2 = 2. / (double)n2 / 9.;
    w3 = 1. - w1;
    w4 = 1. - w2;
    u2 = u * u;
    a = w4 * w4 - u2 * w2;
    b = -2. * w3 * w4;
    c = w3 * w3 - u2 * w1;
    d = b * b - 4 * a * c;
    if(d < 0.)  
        fw = pfsub(a, b, 0.);
    else {
        if(fabs(a) > eps)
            fw = pfsub(a, b, d);
        else {
            if(fabs(b) > eps)   return -c / b;
            fw = pfsub(a, b, 0.);
        }
    }
    return fw * fw * fw;
}
// http://www5.airnet.ne.jp/tomy/cpro/sslib11.htm



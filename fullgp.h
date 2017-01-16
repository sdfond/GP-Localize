const int DATASIZE = 1024;

class FullGP
{
    double getKernelValue(VecDoub loc_a, VecDoub loc_b);
    void addNoise(MatDoub &m);
    void makeKernelVec(VecDoub sample_input[DATASIZE], VecDoub test_point, VecDoub &k_star);
    void makeKernelMat(VecDoub sample_input[DATASIZE], MatDoub &kernel);
    double vecMul(const VecDoub &a, const VecDoub &b);
    VecDoub hyp; //hyper parameters: noise, lenscale1-2, signal, mean

public:
    void init(VecDoub hyper);
    void predict(VecDoub sample[DATASIZE], int sample_size, VecDoub test_input[DATASIZE], VecDoub &pmean, VecDoub &pvar, int test_size);
};

void FullGP::init(VecDoub hyper)
{
    hyp.resize(hyper.size());
    for (int i = 0; i < hyp.size(); i++)
        hyp[i] = hyper[i];
}

void FullGP::predict(VecDoub sample[DATASIZE], int sample_size, VecDoub test_input[DATASIZE], VecDoub &pmean, VecDoub &pvar, int test_size)
{
    MatDoub kernel_mat(0, 0);
    VecDoub y(0), alpha(0), k_star(0), v(0);


    y.resize(sample_size);
    alpha.resize(sample_size);
    k_star.resize(sample_size);
    v.resize(sample_size);
    kernel_mat.resize(sample_size, sample_size);
    for (int i = 0; i < sample_size; i++)
        y[i] = sample[i][2] - hyp[hyp.size()-1];

    makeKernelMat(sample, kernel_mat);
    addNoise(kernel_mat);

    Cholesky ch(kernel_mat);
    ch.solve(y, alpha);

    for (int i = 0; i < test_size; i++)
    {
        makeKernelVec(sample, test_input[i], k_star);
        pmean[i] = vecMul(k_star, alpha) + hyp[hyp.size()-1];
        //predictive mean

        ch.elsolve(k_star, v);
        pvar[i] = getKernelValue(test_input[i], test_input[i]) - vecMul(v, v);
        //predictive variance
    }
}

double FullGP::getKernelValue(VecDoub loc_a, VecDoub loc_b)
{
    double tmp = 0;

    for (int i = 1; i < hyp.size() - 2; i++)
    {
        tmp += pow((loc_a[i-1] - loc_b[i-1]) / hyp[i], 2.0);
    }
    tmp = exp(-0.5 * tmp);

    return hyp[hyp.size()-2] * hyp[hyp.size()-2] * tmp;
}

void FullGP::addNoise(MatDoub &m)
{
    for (int i = 0; i < m.ncols(); i++)
        m[i][i] += hyp[0] * hyp[0];
}

void FullGP::makeKernelVec(VecDoub sample_input[DATASIZE], VecDoub test_point, VecDoub &k_star)
{
    for (int i = 0; i < k_star.size(); i++)
        k_star[i] = getKernelValue(sample_input[i], test_point);
}

void FullGP::makeKernelMat(VecDoub sample_input[DATASIZE], MatDoub &kernel)
{
    for (int i = 0; i < kernel.nrows(); i++)
        for (int j = 0; j < kernel.ncols(); j++)
            kernel[i][j] = getKernelValue(sample_input[i], sample_input[j]);
}

double FullGP::vecMul(const VecDoub &a, const VecDoub &b)
{
    double ret = 0;

    for (int i = 0; i < a.size(); i++)
        ret += a[i] * b[i];

    return ret;
}

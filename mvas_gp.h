#ifndef _MVAS_GP_INCLUDE
#define _MVAS_GP_INCLUDE

#include <cfloat>
#include "nr3/nr3.h"
#include "nr3/cholesky.h"
#include "mvas_util.h"
//#include "mapred_gp.h"
//#undef SEEK_CUR
//#undef SEEK_SET
//#undef SEEK_END
//#include "mpi.h"
#define pic_tset_blk pic_tset_blk4
//@class: covariance funtion
struct mvas_cov
{
    //noise,lenscale #1,...,lenscale #dim_x, signal, mean
    Doub nos;
    VecDoub lsc;
    Doub sig;
    Int dim;
    mvas_cov(VecDoub h, Int d)
    {
        nos = SQR(h[0]);
        lsc.resize(d);
        for (Int i = 1; i <= d; i++)
        {
            lsc[i - 1] = h[i];
        }
        sig = SQR(h[d + 1]);
        dim = d;
    }

    inline Doub se_ard_n(Doub *x, Doub *y)
    {
        //weighted dot product
        Doub val = 0;
        for (Int i = 0; i < dim; i++)
        {
            val += SQR((x[i] - y[i]) / lsc[i]);
            //pcp_s("x:",x[i]);
            //pcp_s("y:",y[i]);
            //pcp_s("l:",lsc[i]);
            //pcp_s("val:",((x[i]-y[i])/lsc[i]));
        }

        //pcp_s("dim:",dim);
        //pcp_s("dist:",val);
        //pcp_s("sig:",sig);
        //pcp_s("nos:",nos);

        //pvec_s("length",lsc);
        return sig * exp(-0.5 * val) + nos;
    }


    inline Doub se_ard(Doub *x, Doub *y)
    {
        //weighted dot product
        Doub val = 0;
        for (Int i = 0; i < dim; i++)
        {
            val += SQR((x[i] - y[i]) / lsc[i]);
            //pcp_s("x:",x[i]);
            //pcp_s("y:",y[i]);
            //pcp_s("l:",lsc[i]);
            //pcp_s("val:",((x[i]-y[i])/lsc[i]));
        }

        //pcp_s("dim:",dim);
        //pcp_s("dist:",val);
        //pcp_s("sig:",sig);
        //pcp_s("nos:",nos);
        //pvec_s("length",lsc);
        return sig * exp(-0.5 * val);
    }

    void se_ard_n(MatDoub a, Int ss, MatDoub &k)
    {
        if (k.nrows() != ss)
        {
            k.resize(ss, ss);
        }
        for (Int i = 0; i < ss; i++)
        {
            for (Int j = i; j < ss; j++)
            {
                //printf("cov:%lf\n",covSEISO(i,j));
                k[i][j] = se_ard(a[i], a[j]);
                if(i == j)
                {
                    k[i][j] += nos;
                }
                else
                {
                    k[j][i] = k[i][j];
                }
            }
        }


    }
    void se_ard(MatDoub a, Int ss, MatDoub &k)
    {

        if (k.nrows() != ss)
        {
            k.resize(ss, ss);
        }

        for (Int i = 0; i < ss; i++)
        {
            for (Int j = i; j < ss; j++)
            {
                //printf("cov:%lf\n",covSEISO(i,j));
                k[i][j] = se_ard(a[i], a[j]);
                if(i != j)
                {
                    k[j][i] = k[i][j];
                }
            }
        }


    }
    //K_AB
    void se_ard(MatDoub a, Int ssa, MatDoub b, Int ssb, MatDoub &k)
    {
        k.resize(ssa, ssb);
        for (Int i = 0; i < ssa; i++)
        {
            for (Int j = 0; j < ssb; j++)
            {

                k[i][j] = se_ard(a[i], b[j]);
            }
        }

    }
};


//@class: Gaussian process and Sparse Gaussian process
struct mvas_gp
{
    mvas_cov *cov;
    Doub h_mu;
    Cholesky *chol_gp; //TODO
    Cholesky *chol_sod;
    //PITC
    Cholesky *chol_pitc_ss;
    Cholesky *chol_pitc_gs;
    MatDoub pitc_aset;
    Int     pitc_anum;

    void init(VecDoub h, Int dx)
    {
        cov = new mvas_cov(h, dx);
        h_mu = h[h.size() - 1];
    }

    inline Cholesky *chol_cov(MatDoub K_dd)
    {
        Cholesky *chol = new Cholesky(K_dd);
        return chol;
    }

    inline Cholesky *chol_cov(MatDoub obs, Int dnum)
    {
        MatDoub K_dd;
        cov->se_ard_n(obs, dnum, K_dd);
        Cholesky *chol = new Cholesky(K_dd);
        return chol;
    }

    inline Cholesky *chol_pcov(MatDoub obs, Int dnum, MatDoub act, Int anum)
    {
        MatDoub kdd;
        Cholesky *chol_kuu = chol_cov(act, anum);
        post_cov(act, anum, chol_kuu, obs, dnum, kdd);
        return chol_cov(kdd);
    }

    inline Cholesky *chol_pcov(MatDoub obs, Int dnum, MatDoub act, Int anum, Cholesky *chol_kuu)
    {
        MatDoub kdd;
        post_cov(act, anum, chol_kuu, obs, dnum, kdd);
        return chol_cov(kdd);
    }

    //@func: d'=L'd where K=LL'; if L is u by u then d is u by n
    inline Int regressor(Cholesky *chol, MatDoub &d)
    {
        Int u = d.nrows();
        for(Int i = 0; i < d.ncols(); i++)
        {
            VecDoub v(u);
            VecDoub alpha(u);
            for(Int j = 0; j < u; j++)
            {
                v[j] = d[j][i];
            }
            chol->solve(v, alpha);
            for(Int j = 0; j < u; j++)
            {
                d[j][i] = alpha[j];
            }
        }

        return SUCC;
    }

    //A: m by n matrix
    //B: n by n matrix
    //c: n-size vector
    //D = A B^{-1} c
    inline Int A_invB_C(MatDoub A, Cholesky *chol_b, VecDoub C, VecDoub &D)
    {
        Int m = A.nrows();
        Int n = A.ncols();
        VecDoub beta(n);
        D.assign(m, 0);
        chol_b->solve(C, beta);
        for(Int k = 0; k < m; k++)
        {
            for(Int j = 0; j < n; j++)
            {
                D[k] += A[k][j] * beta[j];
            }
        }
        return SUCC;

    }


    //A: m by n matrix
    //B: n by n matrix
    //C: n by l matrix
    //D = A B^{-1} C
    inline Int A_invB_C(MatDoub A, Cholesky *chol_b, MatDoub C, MatDoub &D)
    {
        Int m = A.nrows();
        Int n = A.ncols();
        Int l = C.ncols();
        VecDoub v(n);
        VecDoub beta(n);
        VecDoub Cn(n);
        D.assign(m, l, 0);
        for(Int i = 0; i < l; i++)
        {
            //C(:,i)
            for(Int j = 0; j < n; j++)
            {
                Cn[j] = C[j][i];
            }
            //D[i][k]
            chol_b->solve(Cn, beta);
            for(Int k = 0; k < m; k++)
            {
                for(Int j = 0; j < n; j++)
                {
                    D[k][i] += A[k][j] * beta[j];
                }
                //D[i][k] = D[k][i];
            }
        }
        return SUCC;

    }

//A: m by n matrix
    //B: n by n matrix
    //tC: l by n matrix
    //D = A B^{-1} C^T
    inline Int A_invB_Ct(MatDoub A, Cholesky *chol_b, MatDoub C, MatDoub &D)
    {
        Int m = A.nrows();
        Int n = A.ncols();
        Int l = C.nrows();
        VecDoub v(n);
        VecDoub beta(n);
        VecDoub Cn(n);
        D.assign(m, l, 0);
        for(Int i = 0; i < l; i++)
        {
            //C(:,i)
            for(Int j = 0; j < n; j++)
            {
                Cn[j] = C[i][j];
            }
            //D[i][k]
            chol_b->solve(Cn, beta);
            for(Int k = 0; k < m; k++)
            {
                for(Int j = 0; j < n; j++)
                {
                    D[k][i] += A[k][j] * beta[j];
                }
                //D[i][k] = D[k][i];
            }
        }
        return SUCC;

    }


    //A: m by n matrix
    //B: n by n matrix
    //D = A B^{-1} A^T
    inline Int A_invB_C(MatDoub A, Cholesky *chol_b, MatDoub &D)
    {

        Int m = A.nrows();
        Int n = A.ncols();

        VecDoub v(n);
        VecDoub beta(n);
        VecDoub Cn(n);
        D.assign(m, m, 0);

        for(Int i = 0; i < m; i++)
        {

            for(Int j = 0; j < n; j++)
            {
                Cn[j] = A[i][j];
            }
            chol_b->elsolve(Cn, v);
            //D[i][i]
            for(Int j = 0; j < n; j++)
            {
                D[i][i] += v[j] * v[j];
            }
            //D[i][j] i!=j
            chol_b->solve(Cn, beta);
            for(Int k = i + 1; k < m; k++)
            {
                for(Int j = 0; j < n; j++)
                {
                    D[k][i] += A[k][j] * beta[j];
                }
                D[i][k] = D[k][i];
            }
        }

        return SUCC;

    }

    //A: m by n matrix
    //B: n by n matrix
    //C: n by m matrix
    //D = trace(A B^{-1} C)
    inline Int trace_A_invB_C(MatDoub A, Cholesky *chol_b, MatDoub C, VecDoub &D)
    {
        Int m = A.nrows();
        Int n = A.ncols();
        Int l = C.ncols();
        VecDoub v(n);
        VecDoub beta(n);
        VecDoub Cn(n);
        D.assign(m, 0);
        for(Int i = 0; i < l; i++)
        {
            //C(:,i)
            for(Int j = 0; j < n; j++)
            {
                Cn[j] = C[j][i];
            }
            //D[i]
            chol_b->solve(Cn, beta);
            for(Int j = 0; j < n; j++)
            {
                D[i] += A[i][j] * beta[j];
            }
        }
        return SUCC;

    }
    //A: m by n matrix
    //B: n by n matrix
    //D = trace(A B^{-1} A^T)
    inline Int trace_A_invB_C(MatDoub A, Cholesky *chol_b, VecDoub &D)
    {

        Int m = A.nrows();
        Int n = A.ncols();

        VecDoub v(n);
        VecDoub beta(n);
        VecDoub Cn(n);
        D.assign(m, 0);

        for(Int i = 0; i < m; i++)
        {

            for(Int j = 0; j < n; j++)
            {
                Cn[j] = A[i][j];
            }
            chol_b->elsolve(Cn, v);
            //D[i][i]
            for(Int j = 0; j < n; j++)
            {
                D[i] += v[j] * v[j];
            }
        }

        return SUCC;

    }
    //A: m by n matrix
    //B: n by n matrix
    //C: m by n matrix
    //D = trace(A B^{-1} C^T)
    inline Int trace_A_invB_Ct(MatDoub A, Cholesky *chol_b, MatDoub C, VecDoub &D)
    {
        Int m = A.nrows();
        Int n = A.ncols();
        Int l = C.nrows();
        if(m != l)
        {
            throw("something wrong with trace(AB^{-1}C^T)\n");
        }
        VecDoub v(n);
        VecDoub beta(n);
        VecDoub Cn(n);
        D.assign(m, 0);
        for(Int i = 0; i < m; i++)
        {
            //C(:,i)
            for(Int j = 0; j < n; j++)
            {
                Cn[j] = C[i][j];
            }
            //D[i]
            chol_b->solve(Cn, beta);
            for(Int j = 0; j < n; j++)
            {
                D[i] += A[i][j] * beta[j];
            }
        }
        return SUCC;

    }


    //full gp regression
    Int full_reg(MatDoub obs, Int dnum, MatDoub xt, Int ts, VecDoub &t_mu, VecDoub &t_var)
    {
        //step 1:cholesky decomposition
        Int ss = dnum;
        VecDoub v(ss);
        VecDoub alpha(ss);
        //VecDoub beta(ss);
        MatDoub K_dd;
        //build K+\sigma^2_n I.
        cov->se_ard_n(obs, ss, K_dd);
        Cholesky *chol = chol_cov(K_dd);

        //step 2: compute alpha
        for(Int i = 0; i < ss; i++)
        {
            v[i] = obs[i][cov->dim] - h_mu;
        }
        chol->solve(v, alpha);
        //step 3: predictive mean
        //Int ts = xt.nrows();
        t_mu.resize(ts);
        t_var.resize(ts);
        MatDoub K_td;
        cov->se_ard(xt, ts, obs, ss, K_td);

        // \bar{f}_*
        for(Int i = 0; i < ts; i++)
        {
            t_mu[i] = h_mu;
            for(Int j = 0; j < ss; j++)
            {
                t_mu[i] += K_td[i][j] * alpha[j];
            }

            //Step:4 predictive variance
            //k_{i,i}
            t_var[i] = cov->se_ard_n(xt[i], xt[i]);
            //k_{*f}
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = K_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
            chol->elsolve(K_ti, v);
            //predictive variance
            for(Int j = 0; j < ss; j++)
            {
                t_var[i] -= v[j] * v[j];
            }
        }
        delete chol;
        return SUCC;
    }

    // @func  - compute posterior distribution using FGP
    Int full_dist( MatDoub obs, Int ss, Cholesky *chol,
                   MatDoub xt, Int ts,
                   VecDoub &t_mu, MatDoub &t_cov)
    {
        //step 1:cholesky decomposition
        VecDoub v(ss);
        VecDoub vv(ss);
        VecDoub alpha(ss);
        VecDoub beta(ss);
        t_cov.resize(ts, ts);
        t_mu.resize(ts);
        //MatDoub K_dd;
        //build K+\sigma^2_n I.
        //cov->se_ard_n(obs,K_dd);
        //Cholesky chol(K_dd);

        //step 2: compute alpha

        //step 3: predictive mean
        MatDoub K_td;
        cov->se_ard(xt, ts, obs, ss, K_td);
        //pmat_s("ktd",K_td);
        //k_{T,T}
        cov->se_ard_n(xt, ts, t_cov);
        //pmat_s("kss",t_cov);

        for(Int i = 0; i < ss; i++)
        {
            vv[i] = obs[i][cov->dim] - h_mu;
        }
        chol->solve(vv, alpha);

        // \bar{f}_*
        for(Int i = 0; i < ts; i++)
        {

            t_mu[i] = h_mu;
            for(Int j = 0; j < ss; j++)
            {
                t_mu[i] += K_td[i][j] * alpha[j];
            }
            //Step:4 predictive cov
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = K_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
            chol->elsolve(K_ti, v);
            //predictive variance
            for(Int j = 0; j < ss; j++)
            {
                t_cov[i][i] -= v[j] * v[j];
            }
            //var[\hat{f}_ji]
            chol->solve(K_ti, beta);
            for(Int t = i + 1; t < ts; t++)
            {
                for(Int j = 0; j < ss; j++)
                {
                    t_cov[t][i] -= K_td[t][j] * beta[j];
                }
                t_cov[i][t] = t_cov[t][i];
            }
        }
        return SUCC;
    }

    // @func  - compute posterior matrix using FGP
    Int post_var( MatDoub obs, Int ss, Cholesky *chol,
                  MatDoub xt, Int ts,
                  VecDoub &t_var)
    {

        //step 1:cholesky decomposition
        VecDoub v(ss);
        VecDoub beta(ss);
        //MatDoub K_dd;
        //build K+\sigma^2_n I.
        //cov->se_ard_n(obs,K_dd);
        //Cholesky chol(K_dd);

        //step 2: compute alpha

        //step 3: predictive mean
        if (t_var.size() != ts)
        {
            t_var.resize(ts);
        }
        MatDoub K_td;
        cov->se_ard(xt, ts, obs, ss, K_td);
        //pmat_s("ktd",K_td);

        //k_{T,T}
        //pmat_s("test",xt);
        //cov->se_ard_n(xt, ts, t_cov);


        // \bar{f}_*
        for(Int i = 0; i < ts; i++)
        {
            t_var[i] = cov->nos + cov->sig;

            //Step:4 predictive cov
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = K_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
            chol->elsolve(K_ti, v);
            //predictive variance
            for(Int j = 0; j < ss; j++)
            {
                t_var[i] -= v[j] * v[j];
            }
            //var[\hat{f}_ji]
#if 0
            chol->solve(K_ti, beta);
            for(Int t = i + 1; t < ts; t++)
            {
                for(Int j = 0; j < ss; j++)
                {
                    t_cov[t][i] -= K_td[t][j] * beta[j];
                }
                t_cov[i][t] = t_cov[t][i];
            }
#endif
        }
        return SUCC;

    }
    // @func  - compute posterior matrix using FGP
    Int post_cov( MatDoub obs, Int ss, Cholesky *chol,
                  MatDoub xt, Int ts,
                  MatDoub &t_cov)
    {
        //step 1:cholesky decomposition
        VecDoub v(ss);
        VecDoub beta(ss);
        //MatDoub K_dd;
        //build K+\sigma^2_n I.
        //cov->se_ard_n(obs,K_dd);
        //Cholesky chol(K_dd);

        //step 2: compute alpha

        //step 3: predictive mean
        t_cov.resize(ts, ts);
        MatDoub K_td;
        cov->se_ard(xt, ts, obs, ss, K_td);
        //pmat_s("ktd",K_td);

        //k_{T,T}
        //pmat_s("test",xt);
        cov->se_ard_n(xt, ts, t_cov);


        // \bar{f}_*
        for(Int i = 0; i < ts; i++)
        {

            //Step:4 predictive cov
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = K_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
            chol->elsolve(K_ti, v);
            //predictive variance
            for(Int j = 0; j < ss; j++)
            {
                t_cov[i][i] -= v[j] * v[j];
            }
            //var[\hat{f}_ji]
            chol->solve(K_ti, beta);
            for(Int t = i + 1; t < ts; t++)
            {
                for(Int j = 0; j < ss; j++)
                {
                    t_cov[t][i] -= K_td[t][j] * beta[j];
                }
                t_cov[i][t] = t_cov[t][i];
            }
        }
        return SUCC;
    }

    //prepare the local summary for PITC block
    Int pitc_prep( MatDoub D, Int ds, Cholesky *chol_sdd,
                   MatDoub U, Int us, Cholesky *chol_kuu,
                   VecDoub &fu, MatDoub &suu)
    {

        pitc_aset = U;
        pitc_anum = us;
        chol_pitc_ss = chol_kuu;

        //\Sigma_{DD}|U+\sigma_n^2I
        //MatDoub sdd;
        //pcov(U, us, chol_kuu, D, ds, sdd);
        //Cholesky *chol_sdd = new Cholesky(sdd);
        VecDoub v(ds);
        VecDoub alpha(ds);
        VecDoub beta(ds);

        MatDoub K_ud;
        cov->se_ard(U, us, D, ds, K_ud);
        for(Int i = 0; i < ds; i++)
        {
            v[i] = D[i][cov->dim] - h_mu;
        }

        fu.resize(us);
        suu.resize(us, us);
        chol_sdd->solve(v, alpha);
        for(Int i = 0; i < us; i++)
        {
            //step 3: predictive mean
            fu[i] = 0;
            for(Int j = 0; j < ds; j++)
            {
                fu[i] += K_ud[i][j] * alpha[j];
            }
            //Step:4 predictive cov
            VecDoub K_ui(ds);
            for(Int j = 0; j < ds; j++)
            {
                K_ui[j] = K_ud[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
            chol_sdd->elsolve(K_ui, v);
            suu[i][i] = 0;
            for(Int j = 0; j < ds; j++)
            {
                suu[i][i] += v[j] * v[j];
            }
            //var[\hat{f}_ij] where i!=j
            chol_sdd->solve(K_ui, beta);
            for(Int t = i + 1; t < us; t++)
            {
                suu[t][i] = 0;
                for(Int j = 0; j < ds; j++)
                {
                    suu[t][i] += K_ud[t][j] * beta[j];
                }
                suu[i][t] = suu[t][i];
            }
        }
        return SUCC;
    }

    Int pitc_pcov(MatDoub xt, Int ts,  MatDoub &t_cov)
    {

        return pitc_pcov(pitc_aset, pitc_anum,
                         chol_pitc_ss, chol_pitc_gs, xt, ts, t_cov);

    }

    // @func - compute the posterior matrix using PITC
    // @in   - (u,ss) support set
    //       - chol_kuu Cholesky decomp of support set
    //       - chol_suu Cholesky decomp of global summary
    Int pitc_pcov( MatDoub u, Int ss, Cholesky *chol_kuu, Cholesky *chol_suu,
                   MatDoub xt, Int ts,  MatDoub &t_cov)
    {
        //\Sigma_TT
        VecDoub v(ss);
        VecDoub beta(ss);
        post_cov(u, ss, chol_kuu, xt, ts, t_cov);

        MatDoub K_td;
        cov->se_ard(xt, ts, u, ss, K_td);

        for(Int i = 0; i < ts; i++)
        {

            //Step:4 predictive variance
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = K_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
            chol_suu->elsolve(K_ti, v);
            //predictive variance
            for(Int j = 0; j < ss; j++)
            {
                t_cov[i][i] += v[j] * v[j];
            }
            //var[\hat{f}_ij] where j!=i
            chol_suu->solve(K_ti, beta);
            for(Int t = i + 1; t < ts; t++)
            {
                for(Int j = 0; j < ss; j++)
                {
                    t_cov[t][i] += K_td[t][j] * beta[j];
                }
                t_cov[i][t] = t_cov[t][i];
            }
        }

        return SUCC;
    }

    /* specially edited function of post_cov */
    Int post_cov_spec( MatDoub obs, Int ss, Cholesky *chol,
                  MatDoub xt, Int xts, MatDoub yt, Int yts,
                  MatDoub &t_cov)
    {
        //step 1:cholesky decomposition
        VecDoub v(ss);
        VecDoub beta(ss);
        //MatDoub K_dd;
        //build K+\sigma^2_n I.
        //cov->se_ard_n(obs,K_dd);
        //Cholesky chol(K_dd);

        //step 2: compute alpha

        //step 3: predictive mean
        MatDoub Kx_td;
        MatDoub Ky_td;
        cov->se_ard(xt, xts, obs, ss, Kx_td);
        cov->se_ard(yt, yts, obs, ss, Ky_td);
        //pmat_s("ktd",K_td);

        //k_{T,T}
        //pmat_s("test",xt);
        cov->se_ard(xt, xts, yt, yts, t_cov);     //t_cov is resized to xts * yts


        // \bar{f}_*
        for(Int i = 0; i < xts; i++)
        {

            //Step:4 predictive cov
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = Kx_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
//            chol->elsolve(K_ti, v);
            //predictive variance
//            for(Int j = 0; j < ss; j++)
//            {
//                t_cov[i][i] -= v[j] * v[j];
//            }
            //var[\hat{f}_ji]
            chol->solve(K_ti, beta);
            for(Int t = 0; t < yts; t++)
            {
                for(Int j = 0; j < ss; j++)
                {
                    t_cov[i][t] -= Ky_td[t][j] * beta[j];
                }
            }
        }
        return SUCC;
    }

    /* specially edited function of pitc_pcov
       usage: Given data X and data Y, pitc_pcov_spec calculates the predictive covariance between X and Y
       output: t_cov is a matrix of xts * yts
    */
    Int pitc_pcov_spec( MatDoub u, Int ss, Cholesky *chol_kuu, Cholesky *chol_suu,
                        MatDoub xt, Int xts, MatDoub yt, Int yts, MatDoub &t_cov)
    {
        //\Sigma_TT
        VecDoub v(ss);
        VecDoub beta(ss);
        post_cov_spec(u, ss, chol_kuu, xt, xts, yt, yts, t_cov);

        MatDoub Kx_td;
        MatDoub Ky_td;
        cov->se_ard(xt, xts, u, ss, Kx_td);
        cov->se_ard(yt, yts, u, ss, Ky_td);

        for(Int i = 0; i < xts; i++)
        {

            //Step:4 predictive variance
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = Kx_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
//      chol_suu->elsolve(K_ti, v);
//      //predictive variance
//      for(Int j = 0; j < ss; j++) {
//        t_cov[i][i] += v[j] * v[j];
//      }
            //var[\hat{f}_ij] where j!=i
            chol_suu->solve(K_ti, beta);
            for(Int t = 0; t < yts; t++)
            {
                for(Int j = 0; j < ss; j++)
                {
                    t_cov[i][t] += Ky_td[t][j] * beta[j];
                }
            }
        }

        return SUCC;
    }

    // @func - compute the posterior distribution using PITC
    // @in   - (u,ss) support set
    //       - chol_kuu Cholesky decomp of support set
    //       - chol_suu Cholesky decomp of global summary
    Int pitc_dist( MatDoub u, Int ss, Cholesky *chol_kuu,
                   VecDoub g_mu, Cholesky *chol_suu,
                   MatDoub xt, Int ts,
                   VecDoub &t_mu, MatDoub &t_cov)
    {
        //\Sigma_TT
        VecDoub v(ss);
        VecDoub beta(ss);
        VecDoub alpha(ss);
        post_cov(u, ss, chol_kuu, xt, ts, t_cov);
        MatDoub K_td;
        cov->se_ard(xt, ts, u, ss, K_td);
        chol_suu->solve(g_mu, alpha);
        t_mu.resize(ts);
        for(Int i = 0; i < ts; i++)
        {

            // predictive mean for test set
            t_mu[i] = h_mu;
            for(Int j = 0; j < ss; j++)
            {
                t_mu[i] += K_td[i][j] * alpha[j];
            }
            //Step:4 predictive variance
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = K_td[i][j];
            }
            //var[f_ii]
            //v[i]=L\k*[i]
            chol_suu->elsolve(K_ti, v);
            //predictive variance
            for(Int j = 0; j < ss; j++)
            {
                t_cov[i][i] += v[j] * v[j];
            }
            //var[\hat{f}_ij] where j!=i
            chol_suu->solve(K_ti, beta);
            for(Int t = i + 1; t < ts; t++)
            {
                for(Int j = 0; j < ss; j++)
                {
                    t_cov[t][i] += K_td[t][j] * beta[j];
                }
                t_cov[i][t] = t_cov[t][i];
            }
        }

        return SUCC;
    }
  /*
    Int pitc_regr_mr( MatDoub D[], VecInt ds,
                      MatDoub aset, Int as,
                      MatDoub xt, Int ts,
                      VecDoub &t_mu, VecDoub &t_var)
    {

        timespec start, end;
        double elapsed;


        //clock_gettime( 0, &start);
        MatDoub kuu;
        cov->se_ard(aset, as, kuu);
        Cholesky *chol_kuu = chol_cov(kuu);

        int size = ds.size();
        map_data_t mdt[size];

        local_value_t lvt[size];
        global_value_t gvt;
        gvt.gs_kuu = new MatDoub( kuu );
        gvt.gs_zu = new VecDoub( as, (Doub)0 );
        for(Int i = 0; i < size; i++ )
        {
            mdt[i].mvas_pointer = this;
            mdt[i].D = &D[i];
            mdt[i].dsv = ds[i];
            mdt[i].ast = &aset;
            mdt[i].as = as;
            mdt[i].chol_kuu = chol_kuu;

            lvt[i].ls_kuu = new MatDoub();
            lvt[i].ls_zu = new VecDoub();

            mdt[i].lvt = &lvt[i];
            mdt[i].gvt = &gvt;
        }

        mrgp_data_t data;

        data.mdt = mdt;
        data.size = size;
        data.pos = 0;

        //printf( "Before Start !\n" );

        mvasgp_start( &data, 8, size );

        for( Int i = 0; i < size; i++ )
        {
            delete lvt[i].ls_kuu;
            delete lvt[i].ls_zu;
        }

        //clock_gettime( 0, &end);
        elapsed = (double (end.tv_sec - start.tv_sec))  + ((double(end.tv_nsec - start.tv_nsec)) / 1000000000);
        printf( "MR Support Time: %f\n", elapsed );

        //clock_gettime( 0, &start);



        //pitc_regr_low(aset, as, *(gvt.gs_kuu), *(gvt.gs_zu), xt, ts, t_mu, t_var);
        pitc_regr_low2(aset, as, *(gvt.gs_kuu), *(gvt.gs_zu), xt, ts, size, t_mu, t_var);


        //clock_gettime( 0, &end);
        elapsed = (double (end.tv_sec - start.tv_sec))  + ((double(end.tv_nsec - start.tv_nsec)) / 1000000000) ;
        printf( "MR regr_low Time: %f\n", elapsed );
        delete gvt.gs_kuu;
        delete gvt.gs_zu;
        return SUCC;

    }
    static void addVec( double * ls, double * gs, int * len, MPI_Datatype* t)
    {
        cout << "len : " << *len << endl;
        for( int i = 0; i < *len; i++ )
        {
            gs[i] += ls[i];
        }
    }
    static void addMtx( double * ls, double * gs, int * len, MPI_Datatype* t)
    {
        for( int i = 0; i < *len; i++ )
        {
            gs[i] += ls[i];
        }
    }
    Int pitc_regr_mpi(MatDoub D[], VecInt ds,
                      MatDoub aset, Int as,
                      MatDoub xt, Int ts,
                      VecDoub &t_mu, VecDoub &t_var )
    {
        MatDoub kuu;
        cov->se_ard(aset, as, kuu);
        Cholesky *chol_kuu = chol_cov(kuu);

        MatDoub ls_kuu;
        VecDoub ls_zu;
        MatDoub gs_kuu(kuu);
        VecDoub gs_zu(as, (Doub)0);
        MPI_Op op1;
        MPI_Op_create( (MPI_User_function*)addVec, 1, &op1);



        double * lu = new double[as];
        double * lk = new double[as * as];

        double * gu = new double[as];
        double * gk = new double[as * as];

        int r;
        MPI_Comm_rank( MPI_COMM_WORLD, &r);
        Cholesky *chol_sdd = chol_pcov(D[r], ds[r], aset, as, chol_kuu);
        pitc_prep(D[r], ds[r], chol_sdd, aset, as, chol_kuu, ls_zu, ls_kuu);


        for( int i = 0; i < as; i++ )
        {
            lu[i] = ls_zu[i];
            gu[i] = 0;
            for( int j = 0; j < as; j++ )
            {
                lk[i * as + j] = ls_kuu[i][j];
                gk[i * as + j] = 0;
            }
        }


        cout << "Finish local regression " << endl;
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Allreduce( lu, gu, as, MPI_DOUBLE, op1, MPI_COMM_WORLD );
        cout << "Finish global summary mean" << endl;
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Allreduce( lk, gk, as * as, MPI_DOUBLE, op1, MPI_COMM_WORLD );
        cout << "Finish global summary variance" << endl;
        MPI_Barrier( MPI_COMM_WORLD );
        for( int i = 0; i < as; i++ )
        {
            gs_zu[i] = gu[i];
            for( int j = 0; j < as; j++ )
            {
                gs_kuu[i][j] += gk[i * as + j];
            }
        }



        delete lu, delete [] lk, delete [] gu, delete [] gk;

        pitc_regr_low_mpi(aset, as, gs_kuu, gs_zu, xt, ts, ds.size(), t_mu, t_var);
        return SUCC;
    }


    Int pitc_regr_low_mpi( MatDoub u, Int ss,
                           MatDoub suu, VecDoub fu,
                           MatDoub xt, Int ts, Int na,
                           VecDoub &t_mu, VecDoub &t_var)
    {

        Int K = na;
        Int bs = floor(ts / K);
        Int d_xz = xt.ncols();

        MatDoub tset_blk;
        int ts_blk;
        VecDoub pmu_blk;
        VecDoub pvar_blk;

        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );

        if( rank < K - 1 )
        {
            ts_blk = bs;
            tset_blk.resize( bs, d_xz );
            pmu_blk.resize( bs );
            pvar_blk.resize( bs );
            for (Int i = 0; i < bs; i++)
            {
                Int pos = rank * bs + i;
                takeSamp(xt[pos], tset_blk[i], d_xz);
            }
        }
        else
        {
            Int rs = ts - bs * (K - 1);
            ts_blk = rs;
            tset_blk.resize( rs, d_xz );
            pmu_blk.resize( rs );
            pvar_blk.resize( rs );
            for (Int i = 0; i < rs; i++)
            {
                Int pos = bs * (K - 1) + i;
                takeSamp(xt[pos], tset_blk[i], d_xz);
            }
            bs = rs;
        }
        //\Sigma_TT
        VecDoub alpha(ss);
        //MatDoub t_cov;
        MatDoub kuu;
        cov->se_ard(u, ss, kuu);
        //pmat(kuu)r
        Cholesky chol_kuu(kuu);
        //the domain can be extremely large
        //Use post_var instead
        //  this matrix may require a large of memory
        Cholesky chol_suu(suu);

        //t_var.resize(ts);
        chol_suu.solve(fu, alpha);

        // tset_blk, ts_blk, pmu_blk, pvar_blk, alpha,
        // chol_suu, chol_kuu
        cout << "local regression " << endl;
        pitc_regr_low2_core( u, ss, tset_blk, ts_blk, pmu_blk, pvar_blk, alpha, chol_suu, chol_kuu);

        //Use Mapreduce multicore

        if (t_mu.size() != ts)
        {
            t_mu.resize(ts);
        }
        if( t_var.size() != ts )
        {
            t_var.resize(ts);
        }

        double * pmk = new double [bs];
        double * pvk = new double [bs];
        for( int i = 0; i < bs; i++ )
        {
            pmk[i] = pmu_blk[i];
            pvk[i] = pvar_blk[i];
        }


        if( rank == 0 )
        {
            MPI_Status status;
            for( int c = 0; c < bs; c++ )
            {
                t_mu[c] = pmu_blk[c];
                t_var[c] = pvar_blk[c];
            }
            int recSize = 0;
            int start = 0;
            for( int i = 1; i < K; i++ )
            {
                if( i < (K - 1) )
                {
                    recSize = ts / K;
                }
                else
                {
                    recSize = ts - ts / K * i;
                    delete [] pmk;
                    delete [] pvk;
                    pmk = new double[recSize];
                    pvk = new double[recSize];
                }
                start = (ts / K) * i;
                MPI_Recv(pmk, recSize, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status );
                MPI_Recv(pvk, recSize, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status );
                for( int c = 0; c < recSize; c++ )
                {
                    t_mu[start + c] = pmk[c];
                    t_var[start + c] = pvk[c];
                }
            }
        }
        else
        {
            cout << "results send !" << endl;
            MPI_Send(pmk, bs, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );
            MPI_Send(pvk, bs, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );
        }

        delete pmk, delete pvk;
        return SUCC;
    }
  */

// @desc: regression based on PITC
// TODO test the alogrithm
    Int pitc_regr( MatDoub D[], VecInt ds,
                   MatDoub aset, Int as,
                   MatDoub xt, Int ts,
                   VecDoub &t_mu, VecDoub &t_var)
    {


        MatDoub kuu;
        cov->se_ard(aset, as, kuu);
        Cholesky *chol_kuu = chol_cov(kuu);

        MatDoub ls_kuu;
        VecDoub ls_zu;
        MatDoub gs_kuu(kuu);
        VecDoub gs_zu(as, (Doub)0);
        for (Int i = 0; i < ds.size(); i++)
        {
            //local summary in block i
            Cholesky *chol_sdd = chol_pcov(D[i], ds[i], aset, as, chol_kuu);
            pitc_prep(D[i], ds[i], chol_sdd, aset, as, chol_kuu, ls_zu, ls_kuu);
            //global summary
            for(Int r = 0; r < as; r++)
            {
                gs_zu[r] += ls_zu[r];
                for(Int c = 0; c < as; c++)
                {
                    gs_kuu[r][c] += ls_kuu[r][c];
                }
            }
        }
//    pitc_regr_low3(aset, as, gs_kuu, gs_zu, xt, ts, ds.size(), t_mu, t_var);
        return SUCC;

    }
    Int pitc_regr_low3( MatDoub u, Int ss,
                        MatDoub suu, VecDoub fu,
                        MatDoub xt, Int ts, Int na,
                        VecDoub &t_mu, VecDoub &t_var)
    {

        Int K = na;
        Int bs = floor(ts / K);
        Int d_xz = xt.ncols();
        MatDoub *tset_blk;
        tset_blk = new MatDoub[K];
        VecInt ts_blk(K);
        VecDoub *pmu_blk;
        pmu_blk = new VecDoub[K];
        VecDoub *pvar_blk;
        pvar_blk = new VecDoub[K];


        //Assign
        for (Int k = 0; k < K - 1; k++)
        {
            tset_blk[k].resize(bs, d_xz);
            ts_blk[k] = bs;
            pmu_blk[k].resize(bs);
            pvar_blk[k].resize(bs);

            for (Int i = 0; i < bs; i++)
            {
                Int pos = k * bs + i;
                takeSamp(xt[pos], tset_blk[k][i], d_xz);
            }
        }

        //The last block may contain slightly more points
        Int rs = ts - bs * (K - 1);
        tset_blk[K - 1].resize(rs, d_xz);
        ts_blk[K - 1] = rs;
        pmu_blk[K - 1].resize(rs);
        pvar_blk[K - 1].resize(rs);
        for (Int i = 0; i < rs; i++)
        {
            Int pos = bs * (K - 1) + i;
            takeSamp(xt[pos], tset_blk[K - 1][i], d_xz);
        }

        //\Sigma_TT
        VecDoub alpha(ss);
        //MatDoub t_cov;
        MatDoub kuu;
        cov->se_ard(u, ss, kuu);
        //pmat(kuu);
        Cholesky chol_kuu(kuu);
        //the domain can be extremely large
        //Use post_var instead
        //  this matrix may require a large of memory
//    if (flag_for_test == 6 && AZGV == 143)
//    {
//        FILE * fp_test;
//
//        fp_test = fopen("suu.txt", "w");
//        for (int i = 0; i < suu.nrows(); i++)
//        {
//            for (int j = 0; j < suu.ncols(); j++)
//            {
//                fprintf(fp_test, "%lf ", suu[i][j]);
//            }
//            fprintf(fp_test, "\n");
//        }
//    }
        Cholesky chol_suu(suu);

        //t_var.resize(ts);
        chol_suu.solve(fu, alpha);

        // tset_blk, ts_blk, pmu_blk, pvar_blk, alpha,
        // chol_suu, chol_kuu
        for( int i = 0; i < na; i++ )
        {
            pitc_regr_low2_core( u, ss, tset_blk[i], ts_blk[i], pmu_blk[i], pvar_blk[i], alpha, chol_suu, chol_kuu);
        }

        //Use Mapreduce multicore


        /*
           for( int i = 0; i < na; i++ ) {
           post_var(u, ss, &chol_kuu, tset_blk[i], ts_blk[i], pvar_blk[i]);
           MatDoub K_td;
           cov->se_ard( tset_blk[i], ts_blk[i], u, ss, K_td);
           for( int k = 0; k < ts_blk[i]; k++ ) {
           pmu_blk[i][k] = h_mu;
           for(Int j = 0; j < ss; j++) {
           pmu_blk[i][k] += K_td[k][j] * alpha[j];
           }
           VecDoub K_ti(ss);
           for(Int j = 0; j < ss; j++) {
           K_ti[j] = K_td[k][j];
           }
           VecDoub v(ss);
           chol_suu.elsolve(K_ti, v);
           for(Int j = 0; j < ss; j++) {
           pvar_blk[i][k] += v[j] * v[j];
           }
           }
           }*/

        if (t_mu.size() != ts)
        {
            t_mu.resize(ts);
        }
        if( t_var.size() != ts )
        {
            t_var.resize(ts);
        }
        for( int i = 0; i < na; i++ )
        {
            for( int j = 0; j < ts_blk[i]; j++ )
            {
                int ind = i * bs + j;
                t_mu[ind] = pmu_blk[i][j];
                t_var[ind] = pvar_blk[i][j];
            }
        }
        delete [] pmu_blk;
        delete [] pvar_blk;
        delete [] tset_blk;

        return SUCC;
    }

    void pitc_regr_low2_core( MatDoub& u, Int ss, MatDoub & tset_blki, int ts_blki, VecDoub& pmu_blki, VecDoub& pvar_blki, VecDoub& alpha, Cholesky& chol_suu, Cholesky & chol_kuu )
    {
        post_var(u, ss, &chol_kuu, tset_blki, ts_blki, pvar_blki);
        MatDoub K_td;
        cov->se_ard( tset_blki, ts_blki, u, ss, K_td);
        for( int k = 0; k < ts_blki; k++ )
        {
            pmu_blki[k] = h_mu;
            for(Int j = 0; j < ss; j++)
            {
                pmu_blki[k] += K_td[k][j] * alpha[j];
            }
            VecDoub K_ti(ss);
            for(Int j = 0; j < ss; j++)
            {
                K_ti[j] = K_td[k][j];
            }
            VecDoub v(ss);
            chol_suu.elsolve(K_ti, v);
            for(Int j = 0; j < ss; j++)
            {
                pvar_blki[k] += v[j] * v[j];
            }
        }
    }
    //#################end of the structure###################
};
#endif

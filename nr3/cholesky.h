#ifndef NR3_CHOL_INCLUDE
#define NR3_CHOL_INCLUDE
struct Cholesky{
	Int n;
	MatDoub el;
	MatDoub elt;
	Cholesky(MatDoub_I &a) : n(a.nrows()), el(a), elt(a) {
		Int i,j,k;
		VecDoub tmp;
		Doub sum;
		if (el.ncols() != n) throw("need square matrix");
		for (i=0;i<n;i++) {
			for (j=i;j<n;j++) {
				for (sum=el[i][j],k=i-1;k>=0;k--) sum -= el[i][k]*el[j][k];
				if (i == j) {
					if (sum <= 0.0)
						throw("Cholesky failed");
					el[i][i]=sqrt(sum);
				} else el[j][i]=sum/el[i][i];
			}
		}
		for (i=0;i<n;i++) for (j=0;j<i;j++) el[j][i] = 0.;
		for (i=0;i<n;i++) for (j=0;j<n;j++) elt[j][i] = el[i][j];
	}
	void solve(VecDoub_I &b, VecDoub_O &x) {
		Int i,k;
		Doub sum;
		if (b.size() != n || x.size() != n) throw("bad lengths in Cholesky");
		for (i=0;i<n;i++) {
			for (sum=b[i],k=i-1;k>=0;k--) sum -= el[i][k]*x[k];
			x[i]=sum/el[i][i];
		}
		for (i=n-1;i>=0;i--) {
			for (sum=x[i],k=i+1;k<n;k++) sum -= elt[i][k]*x[k];
			x[i]=sum/el[i][i];
		}		
	}
	void elmult(VecDoub_I &y, VecDoub_O &b) {
		Int i,j;
		if (b.size() != n || y.size() != n) throw("bad lengths");
		for (i=0;i<n;i++) {
			b[i] = 0.;
			for (j=0;j<=i;j++) b[i] += el[i][j]*y[j];
		}
	}
	void elsolve(VecDoub_I &b, VecDoub_O &y) {
		Int i,j;
		Doub sum;
		if (b.size() != n || y.size() != n) throw("bad lengths");
		for (i=0;i<n;i++) {
			for (sum=b[i],j=0; j<i; j++) sum -= el[i][j]*y[j];
			y[i] = sum/el[i][i];
		}
	}
	void inverse(MatDoub_O &ainv) {
		Int i,j,k;
		Doub sum;
		ainv.resize(n,n);
		for (i=0;i<n;i++) for (j=0;j<=i;j++){
			sum = (i==j? 1. : 0.);
			for (k=i-1;k>=j;k--) sum -= el[i][k]*ainv[j][k];
			ainv[j][i]= sum/el[i][i];
		}
		for (i=n-1;i>=0;i--) for (j=0;j<=i;j++){
			sum = (i<j? 0. : ainv[j][i]);
			for (k=i+1;k<n;k++) sum -= elt[i][k]*ainv[j][k];
			ainv[i][j] = ainv[j][i] = sum/el[i][i];
		}				
	}

  ///Calculate the determinant of square matrix and retrun in log scale
	Doub logdet() {
		Doub sum = 0.;
		for (Int i=0; i<n; i++) sum += log(el[i][i]);
		return 2.*sum;
	}
};
#endif

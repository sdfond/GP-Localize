struct Expondev : Ran {
	Doub beta;
	Expondev(Doub bbeta, Ullong i) : Ran(i), beta(bbeta) {}
	Doub dev() {
		Doub u;
		do u = doub(); while (u == 0.);
		return -log(u)/beta;
	}
};
struct Logisticdev : Ran {
	Doub mu,sig;
	Logisticdev(Doub mmu, Doub ssig, Ullong i) : Ran(i), mu(mmu), sig(ssig) {}
	Doub dev() {
		Doub u;
		do u = doub(); while (u*(1.-u) == 0.);
		return mu + 0.551328895421792050*sig*log(u/(1.-u));
	}
};
struct Normaldev_BM : Ran {
	Doub mu,sig;
	Doub storedval;
	Normaldev_BM(Doub mmu, Doub ssig, Ullong i)
	: Ran(i), mu(mmu), sig(ssig), storedval(0.) {}
	Doub dev() {
		Doub v1,v2,rsq,fac;
		if (storedval == 0.) {
			do {
				v1=2.0*doub()-1.0;
				v2=2.0*doub()-1.0;
				rsq=v1*v1+v2*v2;
			} while (rsq >= 1.0 || rsq == 0.0);
			fac=sqrt(-2.0*log(rsq)/rsq);
			storedval = v1*fac;
			return mu + sig*v2*fac;
		} else {
			fac = storedval;
			storedval = 0.;
			return mu + sig*fac;
		}
	}
};
struct Cauchydev : Ran {
	Doub mu,sig;
	Cauchydev(Doub mmu, Doub ssig, Ullong i) : Ran(i), mu(mmu), sig(ssig) {}
	Doub dev() {
		Doub v1,v2;
		do {
			v1=2.0*doub()-1.0;
			v2=doub();
		} while (SQR(v1)+SQR(v2) >= 1. || v2 == 0.);
		return mu + sig*v1/v2;
	}
};
struct Normaldev : Ran {
	Doub mu,sig;
	Normaldev(Doub mmu, Doub ssig, Ullong i)
	: Ran(i), mu(mmu), sig(ssig){}
	Doub dev() {
		Doub u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = abs(v) + 0.386595;
			q = SQR(x) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597
			&& (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
		return mu + sig*v/u;
	}
};
struct Gammadev : Normaldev {
	Doub alph, oalph, bet;
	Doub a1,a2;
	Gammadev(Doub aalph, Doub bbet, Ullong i)
	: Normaldev(0.,1.,i), alph(aalph), oalph(aalph), bet(bbet) {
		if (alph <= 0.) throw("bad alph in Gammadev");
		if (alph < 1.) alph += 1.;
		a1 = alph-1./3.;
		a2 = 1./sqrt(9.*a1);
	}
	Doub dev() {
		Doub u,v,x;
		do {
			do {
				x = Normaldev::dev();
				v = 1. + a2*x;
			} while (v <= 0.);
			v = v*v*v;
			u = doub();
		} while (u > 1. - 0.331*SQR(SQR(x)) &&
			log(u) > 0.5*SQR(x) + a1*(1.-v+log(v)));
		if (alph == oalph) return a1*v/bet;
		else {
			do u=doub(); while (u == 0.);
			return pow(u,1./oalph)*a1*v/bet;
		}
	}
};

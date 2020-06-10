#include "iml.h"

bool iml::BICGSTAB(const Qcx_mat & _A, const Qcx_vec & _b, Qcx_vec & _x,
                        value_t _eps, int _max_iter, int& _iterNum)
{
    _iterNum = 0;
    _x.zeros(arma::size(_b));
    Qcx_vec r0 = _b - _A * _x;
    Qcx_vec r1 = r0;
    Complex rou0(1, 0), alph(1, 0), omiga(1, 0);
    Qcx_vec v, p0, p1, s, t;
    v.zeros(arma::size(_x));
    p0.zeros(arma::size(_x));
    Complex rou1, beta;
    auto b_normal = arma::norm(_b);
    for (int i = 0; i < _max_iter; ++i)
    {
        ++_iterNum;
        rou1 = arma::cdot(r1, r0);
        beta = (rou1 / rou0) * (alph / omiga);
        p1 = r0 + beta * (p0 - omiga * v);
        v = _A * p1;
        alph = rou1 / (arma::cdot(r1, v));
        _x = _x + alph * p1;
        s = r0 - alph * v;
        auto resdual = arma::norm(s) / b_normal;
		Qcout << std::setw(15) << i << std::setw(15) << resdual << std::endl;
        if (resdual < _eps) break;
        t = _A * s;
        omiga = arma::cdot(t, s) / arma::cdot(t, t);
        _x = _x + omiga * s;
        r0 = s - omiga * t;
        resdual = arma::norm(s) / b_normal;
        if (resdual < _eps) break;
        rou0 = rou1;
        p0   = p1;
    }
    return _iterNum != _max_iter;
}

bool iml::BICGSTAB(const CMatrix& _A, const CVector& _b, CVector& _x, value_t _eps, int _max_iter, int& _iterNum)
{
	_iterNum = 0;
	//_x.zeros(arma::size(_b));
	_x.setZero(_b.size());
	//Qcx_vec r0 = _b - _A * _x;
	CVector r0= _b - _A * _x;
	CVector r1 = r0;
	Complex rou0(1, 0), alph(1, 0), omiga(1, 0);
	CVector v, p0, p1, s, t;
	//v.zeros(arma::size(_x));
	//p0.zeros(arma::size(_x));
	v.setZero(_x.size());
	p0.setZero(_x.size());
	Complex rou1, beta;
	//auto b_normal = arma::norm(_b);
	auto b_normal = _b.norm();
	for (int i = 0; i < _max_iter; ++i)
	{
		++_iterNum;
		//rou1 = arma::cdot(r1, r0);
		rou1 = r1.dot(r0);
		beta = (rou1 / rou0) * (alph / omiga);
		p1 = r0 + beta * (p0 - omiga * v);
		v = _A * p1;
		//alph = rou1 / (arma::cdot(r1, v));
		alph = rou1 / (r1.dot(v));
		_x = _x + alph * p1;
		s = r0 - alph * v;
		//auto resdual = arma::norm(s) / b_normal;
		auto resdual = s.norm() / b_normal;
		Qcout << std::setw(15) << i << std::setw(15) << resdual << std::endl;
		if (resdual < _eps) break;
		t = _A * s;
		//omiga = arma::cdot(t, s) / arma::cdot(t, t);
		omiga = t.dot(s) / t.dot(t);
		_x = _x + omiga * s;
		r0 = s - omiga * t;
		//resdual = arma::norm(s) / b_normal;
		resdual = s.norm() / b_normal;
		if (resdual < _eps) break;
		rou0 = rou1;
		p0 = p1;
	}
	return _iterNum != _max_iter;
}

bool iml::BICGSTAB(CVector& x, MulVct mul, const CVector & b,int iterMax, value_t iterEps, int &iterNum)
{
	auto n = b.size();
	x = CVector::Zero(n);

	CVector r = b - mul(x);
	CVector r0 = r;
	value_t r0_sqnorm = r0.squaredNorm();
	value_t rhs_sqnorm = b.squaredNorm();

	if (rhs_sqnorm == 0)
		return true;

	Complex rho = 1;

	Complex alpha = 1;
	Complex w = 1;

	CVector v = CVector::Zero(n), p = CVector::Zero(n);
	CVector s(n), t(n);

	value_t tol2 = iterEps * iterEps;
	value_t eps2 = Eigen::NumTraits<Complex>::epsilon()*
		Eigen::NumTraits<Complex>::epsilon();
	int i = 0;
	int restarts = 0;

	while (r.squaredNorm() / rhs_sqnorm > tol2 && i < iterMax)
	{
		Complex rho_old = rho;

		rho = r0.dot(r);
		if (abs(rho) < eps2*r0_sqnorm)
		{
			r0 = r;
			rho = r0_sqnorm = r.squaredNorm();
			if (restarts++ == 0)
				i = 0;
		}
		Complex beta = (rho / rho_old) * (alpha / w);
		p = r + beta * (p - w * v);

		v.noalias() = mul(p);

		alpha = rho / r0.dot(v);
		s = r - alpha * v;

		t.noalias() = mul(s);

		value_t tmp = t.squaredNorm();

		if (tmp > value_t(0))
			w = t.dot(s) / tmp;
		else
			w = value_t(0);
		x += alpha * p + w * s;
		r = s - w * t;

		Qcout << ++i << " " << sqrt(r.squaredNorm() / rhs_sqnorm) << std::endl;
	}
	value_t tol_error = (r.squaredNorm() / rhs_sqnorm);
	//if (tol_error > tol2)
		//return false;

	Qcout << "Iteration Num:" << i << std::endl;
	return true;
}

bool iml::GMRES(const CMatrix& _A, const CVector& _b, CVector& _x, value_t _eps, int _max_iter,int _out_iter)
{
	CMatrix v,h;
	CVector e,yk,x0,xk,r0;
	int m;

	v.resize(_A.rows(), _max_iter + 1);
	h.setZero(_max_iter + 1, _max_iter);
	e.setZero(_max_iter+1);
	x0.setZero(_A.rows());
	int cur = 0;

	
	while (true)
	{
		cur++;
		m = _max_iter;
		r0 = _b - _A * x0;
		auto beta = r0.norm();
		v.col(0) = r0 *(1.0f/ beta);

		Qcout << std::setw(15) << cur << std::setw(15) << r0.norm() << std::endl;
		if (r0.norm() < _eps)
		{
			_x = x0;
			break;
		}
		
		for (int j = 0; j < _max_iter; j++)
		{
			CVector w = _A * v.col(j);
			for (int i = 0; i < j+1; i++)
			{
				h(i, j) = w.adjoint()*v.col(i);
				w = w - (h(i, j)*v.col(i));
			}
			//CVector v_temp = Av - v.leftCols(j)*(h.col(j).segment(0, j));
			h(j + 1, j) = w.norm();
			if (w.norm() < 1e-6)
			{
				m = j+1;
				break;
			}
			v.col(j + 1) = w *(1.0f/ h(j + 1, j));
		}
		e(0) = beta;
		yk = h.block(0,0,m+1,m).colPivHouseholderQr().solve(e);
		x0 = x0 + v.leftCols(m) * yk;
	}
	return cur != _out_iter;
	
}


bool iml::CG(const Qcx_mat & _A, const Qcx_vec & _b, Qcx_vec & _x, value_t _eps, int _max_iter, int & _iterNum)
{
    _iterNum = 0;
    _x.zeros(_A.n_rows);
	Qcx_vec r = _b - _A * _x;
    auto p = r;
    auto rsold = arma::dot(r, r);
    while (_iterNum < _max_iter)
    {
        ++_iterNum;
		Qcx_vec Ap = _A * p;
        auto alpha = rsold / (arma::dot(p, Ap));
        _x = _x + alpha * p;
        r = r - alpha * Ap;
        auto rsnew = arma::dot(r, r);
		Qcout << std::setw(15) << _iterNum << std::setw(15) << std::sqrt(std::abs(rsnew)) << std::endl;
        if (std::sqrt(std::abs(rsnew)) < _eps) break;
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    }
    return _iterNum != _max_iter;
}

bool iml::CG(const CMatrix& _A, const CVector& _b, CVector& _x, value_t _eps, int _max_iter, int& _iterNum)
{
	_iterNum = 0;
	_x.setZero(_A.rows());
	CVector r0 = _b - _A * _x;
	CVector p = r0;
	auto rsold = r0.dot(r0);
	while (_iterNum < _max_iter)
	{
		++_iterNum;
		CVector Ap = _A * p;
		Complex alpha = ((r0.transpose()*r0) / (p.transpose()*Ap))(0,0);
		//Qcout << r0.dot(r0) <<'    '<< p.dot(Ap)<< std::endl;
		//system("pause");
		_x = _x + alpha * p;
		CVector r1 = r0 - alpha * Ap;
		Qcout << std::setw(15) << _iterNum << std::setw(15) << r1.norm() << std::endl;
		if (r1.norm() < _eps)
			break;
		auto beta = ((r1.transpose()*r1) / (r0.transpose()*r0))(0,0);
		p = r1 + beta * p;
		r0 = r1;
		
	}
	return _iterNum != _max_iter;
}

bool iml::BICG(const Qcx_mat & _A, const Qcx_vec & _b, Qcx_vec & _x, value_t _eps, int _max_iter, int & _iterNum)
{
    return false;
}

bool iml::CGNE(const Qcx_mat & _A, const Qcx_vec & _b, Qcx_vec & _x, value_t _eps, int _max_iter, int & _iterNum)
{
    _iterNum = 0;
    _x.zeros(_A.n_cols);
    Qcx_vec r0 = _b - _A * _x;
    Qcx_vec p0 = (r0.t() * _A).t();
    Qcx_vec p1;
    while (_iterNum < _max_iter)
    {
        ++_iterNum;
        auto alpha = (arma::cdot(r0, r0) / arma::cdot(p0, p0)).real();
        _x = _x + alpha * p0;
        r0 = r0 - alpha * (_A * p0);
        auto rs = arma::norm(r0);
        if (rs < _eps) break;
        auto pn2 = arma::cdot(p0, p0);
        p1 = (r0.t() * _A).st();
        auto beta = -(arma::dot(p1, p0) / pn2);
        p0 = arma::conj(p1) + beta * p0;
    }
    return _iterNum != _max_iter;
}

bool iml::GCR(const Qcx_mat & A, const Qcx_vec & b, Qcx_vec & x, value_t threshold, int max_iter, int& iter_num)
{
    x.zeros(b.size());
    value_t b_norm = arma::norm(b);
    Qcx_vec r0 = b - A * x;
    Qcx_vec p0 = r0;
    Qcx_vec ar0 = A * r0;
    Qcx_vec ap0 = A * p0;

    while (iter_num < max_iter)
    {
        ++iter_num;
        auto r0ar0 = arma::cdot(r0, ar0);
        auto alpha = r0ar0 / arma::cdot(ap0, ap0);
        x = x + alpha * p0;
        r0 = r0 - alpha * ap0;
        auto resdual = arma::norm(r0) / b_norm;
        if (resdual < threshold) break;
        ar0 = A * r0;
        auto beta = arma::cdot(r0, ar0) / r0ar0;
        p0 = r0 + beta * p0;
        ap0 = ar0 + beta * ap0;
    }
    return iter_num != max_iter;
}

///////////////////////////////////////////////////////////////////////////////

void test::GCR()
{
    Qcx_mat A = arma::randu<Qcx_mat>(10, 10);
    A = A.t() * A;
    Qcx_vec b = arma::randu<Qcx_vec>(10);
    Qcx_vec x0 = arma::solve(A, b);
    x0.print("exact x:");
    int iter_num = 0;
    iml::GCR(A, b, x0, 3e-3f, 500, iter_num);
    x0.print("GCR x:");
}

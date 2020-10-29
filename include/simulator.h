#ifndef __SIMULATOR__
#define __SIMULATOR__

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <mesh_object.h>

//Simulator Base Class
 class Simulator{
    protected:
        Eigen::SparseMatrix<double> M, M_inverse;
        Eigen::VectorXd q_curr, q_prev, q_dot_curr, q_dot_prev, L;
        Eigen::SparseMatrix<double> P;

        const MeshObject& obj;
        const double dt;
        const double k; //spring constant

        Simulator(const MeshObject& obj, double k, double dt);
        void d2Vdq2(Eigen::SparseMatrix<double>& d2Vdq2, const Eigen::VectorXd& q_curr) const;
        void dVdq(Eigen::VectorXd& dVdq, const Eigen::VectorXd& q_curr) const;
        void f_ext(Eigen::VectorXd& f_ext);

    public:
        virtual ~Simulator();
        virtual void set_initial_config();
        virtual void step()=0;
        Eigen::Map<const Eigen::MatrixXd> V() const; //behaves similarly to a constant reference
};

#endif //__SIMULATOR__
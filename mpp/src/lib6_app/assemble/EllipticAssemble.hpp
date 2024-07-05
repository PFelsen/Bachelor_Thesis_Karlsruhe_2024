#ifndef ELLIPTIC_ASSEMBLE_HPP
#define ELLIPTIC_ASSEMBLE_HPP

#include "EllipticProblems.hpp"
#include "IEllipticAssemble.hpp"
#include "ElementPool.hpp"

template<typename ElementType, typename FaceElementType>
class EllipticAssemble: public IEllipticAssemble {
protected:
    mutable ElementPool<ElementType> elementPool;
    mutable FaceElementPool<FaceElementType> faceElementPool;

public: 
    explicit EllipticAssemble(const IEllipticProblem &problem) : IEllipticAssemble(problem) {};

    double Energy(const Vector &u) const override {
        double energy = 0.0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        const auto &elem = elementPool.Get(u, *c);

        for (int q = 0; q < elem.nQ(); ++q) {
            double w = elem.QWeight(q);
            VectorField DU = elem.Derivative(q, u);
            Tensor K = problem.Permeability(elem.QPoint(q));
            VectorField K_DU = K * DU;
            energy += w * K_DU * DU;
        }
        }
        return sqrt(PPM->SumOnCommSplit(energy, u.CommSplit()));
    }

    double H1(const Vector &u) const override {
        double err = 0.0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        const auto &elem = elementPool.Get(u, *c);

        for (int q = 0; q < elem.nQ(); ++q) {
            double w = elem.QWeight(q);
            Scalar U = elem.Value(q, u);
            VectorField DU = elem.Derivative(q, u);
            err += w * U * U + w * DU * DU;
        }
        }
        return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
    }

    double EnergyError(const Vector &u) const override {
        double err = 0.0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        const auto &elem = elementPool.Get(u, *c);

        for (int q = 0; q < elem.nQ(); ++q) {
            double w = elem.QWeight(q);
            Point z = elem.QPoint(q);
            Tensor K = problem.Permeability(z);
            Tensor IK = Invert(K);
            VectorField diff = (K * elem.Derivative(q, u) - problem.Flux(z));
            err += w * (IK * diff) * diff;
        }
        }
        return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
    }

    double L2(const Vector &u) const override {
        double l2 = 0.0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        const auto &elem = elementPool.Get(u, *c);

        for (int q = 0; q < elem.nQ(); ++q) {
            double w = elem.QWeight(q);
            Scalar U = elem.Value(q, u);
            l2 += w * U * U;
        }
        }
        return sqrt(PPM->SumOnCommSplit(l2, u.CommSplit()));
    }

    double L2Error(const Vector &u) const override {
        double err = 0.0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        const auto &elem = elementPool.Get(u, *c);

        for (int q = 0; q < elem.nQ(); ++q) {
            double w = elem.QWeight(q);
            Scalar U = elem.Value(q, u);
            Scalar Sol = problem.Solution(elem.QPoint(q));
            err += w * (U - Sol) * (U - Sol);
        }
        }
        return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
    }

    double L2CellAvgError(const Vector &u) const override {
        double err = 0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        const auto &elem = elementPool.Get(u, *c);

        double w = 0;
        Scalar U = 0;
        Scalar Sol = 0;
        for (int q = 0; q < elem.nQ(); ++q) {
            w += elem.QWeight(q);
            U += elem.Value(q, u);
            Sol += problem.Solution(elem.QPoint(q));
        }
        err += w * (U - Sol) * (U - Sol);
        }
        return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
    }

    double MaxError(const Vector &u) const override {
        double err = 0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        const auto &elem = elementPool.Get(u, *c);

        for (int q = 0; q < elem.nQ(); ++q) {
            Scalar p = elem.Value(q, u);
            err = max(err, abs(p - problem.Solution(elem.QPoint(q))));
        }
        }
        return PPM->Max(err, u.CommSplit());
    }

    double FaceError(const Vector &u) const override {
        double face_error = 0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        for (int face = 0; face < c.Faces(); ++face) {
            const auto &faceElem = faceElementPool.Get(u, *c, face);

            for (int q = 0; q < faceElem.nQ(); ++q) {
            double w = faceElem.QWeight(q);
            Scalar U = faceElem.Value(q, u);
            Scalar Sol = problem.Solution(faceElem.QPoint(q));
            face_error += w * (U - Sol) * (U - Sol);
            }
        }
        }
        return sqrt(PPM->SumOnCommSplit(face_error, u.CommSplit()));
    }

    double FluxError(const Vector &u) const override {
        double flux_error = 0;
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
        BFParts bnd(u.GetMesh(), *c);
        if (!bnd.onBnd()) continue;
        for (int face = 0; face < c.Faces(); ++face) {
            if (bnd[face] == 2) {
            const auto &faceElem = faceElementPool.Get(u, *c, face);

            for (int q = 0; q < faceElem.nQ(); ++q) {
                double w = faceElem.QWeight(q);
                VectorField G = problem.Flux(faceElem.QPoint(q));
                Tensor K = problem.Permeability(faceElem.QPoint(q));
                VectorField DU = faceElem.Derivative(q, u);
                Scalar F = (K * DU - G) * faceElem.QNormal(q);
                flux_error += w * F * F;
            }
            } else if (bnd[face] == 0) {
            const auto &faceElem = faceElementPool.Get(u, *c, face);

            for (int q = 0; q < faceElem.nQ(); ++q) {
                double w = faceElem.QWeight(q);
                Tensor K = problem.Permeability(faceElem.QPoint(q));
                VectorField DU = faceElem.Derivative(q, u);
                Scalar F = K * DU * faceElem.QNormal(q);
                flux_error += w * F * F;
            }
            }
        }
        }
        return sqrt(PPM->SumOnCommSplit(flux_error, u.CommSplit()));
    }
};

#endif // ELLIPTIC_ASSEMBLE_HPP
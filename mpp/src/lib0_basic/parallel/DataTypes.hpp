#ifndef CARDMECH_DATATYPES_HPP
#define CARDMECH_DATATYPES_HPP

#include <complex>
#include <mpi.h>

template<typename T>
struct has_mpi_datatype :
    std::integral_constant<
        bool, std::is_same_v<char, typename std::remove_cv<T>::type>
                  || std::is_same_v<short, typename std::remove_cv<T>::type>
                  || std::is_same_v<int, typename std::remove_cv<T>::type>
                  || std::is_same_v<long, typename std::remove_cv<T>::type>
                  || std::is_same_v<unsigned long, typename std::remove_cv<T>::type>
                  || std::is_same_v<float, typename std::remove_cv<T>::type>
                  || std::is_same_v<double, typename std::remove_cv<T>::type>
                  || std::is_same_v<std::complex<double>, typename std::remove_cv<T>::type>> {};

template<typename T>
struct mpi_datatype {
  static MPI_Datatype value() noexcept { return MPI_DATATYPE_NULL; }
};

template<>
struct mpi_datatype<char> {
  static MPI_Datatype value() noexcept { return MPI_CHAR; }
};

template<>
struct mpi_datatype<short> {
  static MPI_Datatype value() noexcept { return MPI_SHORT; }
};

template<>
struct mpi_datatype<int> {
  static MPI_Datatype value() noexcept { return MPI_INT; }
};

template<>
struct mpi_datatype<long> {
  static MPI_Datatype value() noexcept { return MPI_LONG; }
};

template<>
struct mpi_datatype<unsigned long> {
  static MPI_Datatype value() noexcept { return MPI_UNSIGNED_LONG; }
};

template<>
struct mpi_datatype<float> {
  static MPI_Datatype value() noexcept { return MPI_FLOAT; }
};

template<>
struct mpi_datatype<double> {
  static MPI_Datatype value() noexcept { return MPI_DOUBLE; }
};

template<>
struct mpi_datatype<std::complex<double>> {
  static MPI_Datatype value() noexcept { return MPI_DOUBLE_COMPLEX; }
};

#endif // CARDMECH_DATATYPES_HPP

// class for objects that represent scalars
template <typename T>
class A_Scalar {
  private:
    T const& s;  // value of the scalar

  public:
    mutable bool iscached = 1;
    bool evalme = 0;
    // constructor initializes value
    A_Scalar (T const& v)
     : s(v) {
    }

    // for index operations the scalar is the value of each element
    T operator[] (size_t) const {
        return s;
    }


    void cache() const {}

    T val() const {
      return s;
    }
    
    // for index operations the scalar is the value of each element
    // match () 1D access.
    T operator() (size_t i)  const {
        return s;
    }
    // match () 2D access..
    T operator() (size_t i, size_t j)  const {
        return s;
    }

    // scalars have zero as size
    size_t size() const {
        return 0;
    };
    size_t getnrows() const {
        return 0;
    };
    size_t getncols() const {
        return 0;
    };
};

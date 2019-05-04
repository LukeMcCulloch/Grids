


template <typename LHS, typename RHS>
class MatrixSum{

  public:
    using value_type = typename LHS::value_type;

    MatrixSum(const LHS& lhs, const RHS& rhs) : rhs(rhs), lhs(lhs) {}

    // compute sum when value requested
    value_type operator() (int x, int y) const  {
        return lhs(x, y) + rhs(x, y);
    }
  private:
    const LHS& lhs;
    const RHS& rhs;
};


template <typename LHS, typename RHS>
MatrixSum<LHS, RHS> operator+(const LHS& lhs, const LHS& rhs) {
    return MatrixSum<LHS, RHS>(lhs, rhs);
}
template<typename I1, typename I2, typename I3,
   typename O1, typename O2>


class BaseTest01{

public:
   BaseTest01(){}
   bool run(int first, int last, std::fstream &fout);
   virtual void
   apply(const I1 &, const I2 &, const I3 &, O1 &, O2&)=0;
   virtual const char *getName() const=0;
protected:
   void
   addCase(const I1 &i1, const I2 &i2, const I3 &i3,
      const O1 &o1, const O2 &o2);
private:
   std::vector<I1> i1s;
   std::vector<I2> i2s;
   std::vector<I3> i3s;
   std::vector<O1> o1s;
   std::vector<O2> o2s;
};


template<typename I1, typename I2, typename I3,
   typename O1, typename O2>
bool
BaseTest01<I1, I2, I3, O1, O2>::run(int first, int last,
   std::fstream &fout)
{
   O1 tempo1;
   O2 tempo2;

   // constrain range, if necessary
   if(first < 0)
      first = 0;
   if(last > i1s.size() - 1)
      last = i1s.size() - 1;

   for(int i = first; i <= last; ++i)
   {
      // call polymorphically
      apply(i1s[i], i2s[i], i3s[i], tempo1, tempo2);
      if((tempo1 == o1s[i]) && (tempo2 == o2s[i]))
         fout << getName()
              << "(" << i << ") successful. Output: "
              << tempo1 << ", " << tempo2 << std::endl;
      else
      {
         fout << getName()
              << "(" << i << ") UNSUCCESSFUL" << std::endl;
         fout << "  expected: " << o1s[i]
              << ", actual: " << tempo1 << std::endl;
         fout << "  expected: " << o2s[i]
              << ", actual: " << tempo2 << std::endl;
         return false;
      }
   }
   return true;
}
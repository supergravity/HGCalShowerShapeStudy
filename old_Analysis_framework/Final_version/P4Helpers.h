#ifndef P4HELPERS
#define P4HELPERS

#include <vector>
#include <limits>

namespace P4Helpers
{
 
  /// Computes efficiently \Delta{\eta}
  //Delta Eta
  template <typename T , typename U>
    double deltaEta( const T& p1, const U& p2 )
    {
      return p1.Eta()-p2.Eta();
    }

  //Delta Phi
  template <typename T , typename U>
    double deltaPhi( const T& pA, const  U& pB )
    {
     return  -remainder( -pA.phi() + pB.phi(), 2*M_PI );
    }

  //Delta R
  template <typename T , typename U>
    double deltaR( const T& pA, const U& pB )
    {
      //using std::abs;
      //using std::sqrt;
      const double dEta = P4Helpers::deltaEta( pA, pB );
      const double dPhi = P4Helpers::deltaPhi( pA, pB );
      return sqrt( dEta*dEta + dPhi*dPhi );
    }


   //Invariant Mass
   template <typename T , typename U>
   double invMass( const T& pA, const U& pB )
   { 
    return (pA+pB).M(); 
   }
  
  template <typename T , typename U , typename A > 
    bool closestDeltaR( const T& pA, std::vector<U> & coll, A& index, double& deltaR )
    {
  
      //deltaR = std::numeric_limits<double>::max(); // big value
      deltaR = std::numeric_limits<double>::max(); // big value
      bool l_return = false;
      A l_idx = 0;
      for (unsigned int i=0; i<coll.size(); ++i,++l_idx) {
	double rtu = P4Helpers::deltaR(pA,coll.at(i));
	if ( rtu < deltaR ) {
	  index  = l_idx;
	  deltaR = rtu;
	  l_return = true;
	}
      }
      return l_return;
    }

}
#endif

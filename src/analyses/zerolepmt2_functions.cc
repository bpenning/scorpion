#include "zerolepmt2_functions.hh"
std::vector<int> getPseudoJetsGrouping(std::vector<jjet> Object){
    std::vector<int> Object_Group; 
    std::vector<double> Axis1;
    std::vector<double> Axis2;
    int vsize = (int) Object.size();

  
  Object_Group.clear();
  Axis1.clear();
  Axis2.clear();

  for(int j = 0; j < vsize; j++){
    Object_Group.push_back(0);
  }

  for(int j = 0; j < 5; j++){
    Axis1.push_back(0);
    Axis2.push_back(0);
  }
  
 
  for (int i = 0; i <vsize; i++){
   
    if ( (Object[i]).P() > (Object[i]).E() + 0.001) { 
 
        std::cout << "Object " << i << " has E = " << (Object[i]).E()
                                        << " less than P = " << (Object[i]).P() ;
    
    } 
  } 

   
   
  int I_Max = -1;
  int J_Max = -1;
   

    double Mass_Max = 0.;
    double InvariantMass = 0.;
    
    // maximize the invariant mass of two objects
    for (int i=0;i<vsize;i++){    
      Object_Group[i] = 0;
        for (int j=i+1;j<vsize;j++){  

            // either the invariant mass
            InvariantMass =  ((Object[i]).E() +  (Object[j]).E())* ((Object[i]).E() + (Object[j]).E())
              - ((Object[i]).Px() + (Object[j]).Px())* ((Object[i]).Px() + (Object[j]).Px()) 
              - ((Object[i]).Py() + (Object[j]).Py())* ((Object[i]).Py() + (Object[j]).Py())
              - ((Object[i]).Pz() + (Object[j]).Pz())* ((Object[i]).Pz() + (Object[j]).Pz()) ;  
            if ( Mass_Max < InvariantMass){
              Mass_Max = InvariantMass;
              I_Max = i;
              J_Max = j;
            }
        }
    }
    
    if (J_Max>0) {

      Axis1[0] = (Object[I_Max]).Px() /  (Object[I_Max]).P();
      Axis1[1] = (Object[I_Max]).Py() /  (Object[I_Max]).P();
      Axis1[2] = (Object[I_Max]).Pz() /  (Object[I_Max]).P();
    
      Axis1[3] = (Object[I_Max]).P();
      Axis1[4] = (Object[I_Max]).E(); 
  
      Axis2[0] = (Object[J_Max]).Px() /  (Object[J_Max]).P();
      Axis2[1] =(Object[J_Max]).Py() /  (Object[J_Max]).P();
      Axis2[2] = (Object[J_Max]).Pz() /  (Object[J_Max]).P();
    
      Axis2[3] = (Object[J_Max]).P();
      Axis2[4] = (Object[J_Max]).E(); 

    } else {
        std::cout << "ERROR: in step 1 of the algorithm; line "<< __LINE__ << std::endl;
    }
    
    
   
  // seeding done 
  // now do the hemisphere association
   

  int numLoop = 0;
  bool I_Move = true;


  while (I_Move && (numLoop < 100)){

    I_Move = false;
    numLoop++;
   
    double Sum1_Px = 0.;
    double Sum1_Py = 0.;
    double Sum1_Pz = 0.;
    double Sum1_P = 0.;
    double Sum1_E = 0.; 
    double Sum2_Px = 0.;
    double Sum2_Py = 0.;
    double Sum2_Pz = 0.;
    double Sum2_P = 0.;
    double Sum2_E = 0.; 
   
    
    
      for (int i=0;i<vsize;i++){  
        if (i == I_Max) {
	  Object_Group[i] = 1;
	  Sum1_Px += (Object[i]).Px();
	  Sum1_Py += (Object[i]).Py();
	  Sum1_Pz += (Object[i]).Pz();
	  Sum1_P += (Object[i]).P();
	  Sum1_E += (Object[i]).E(); 
	} else if (i == J_Max) {
	  Object_Group[i] = 2;
	  Sum2_Px += (Object[i]).Px();
	  Sum2_Py += (Object[i]).Py();
	  Sum2_Pz += (Object[i]).Pz();
	  Sum2_P += (Object[i]).P();
	  Sum2_E += (Object[i]).E(); 
        } else {
	
	
          if(!I_Move){ 
	  
            double NewAxis1_Px = Axis1[0] * Axis1[3];
            double NewAxis1_Py = Axis1[1] * Axis1[3];
            double NewAxis1_Pz = Axis1[2] * Axis1[3];
            double NewAxis1_E = Axis1[4];
	 
            double NewAxis2_Px = Axis2[0] * Axis2[3];
            double NewAxis2_Py = Axis2[1] * Axis2[3];
            double NewAxis2_Pz = Axis2[2] * Axis2[3];
            double NewAxis2_E = Axis2[4];
       
               
	  
            double mass1 =  NewAxis1_E - (((Object[i]).Px()*NewAxis1_Px + (Object[i]).Py()*NewAxis1_Py +
                                          (Object[i]).Pz()*NewAxis1_Pz)/(Object[i]).P());
	 
            double mass2 =  NewAxis2_E - (((Object[i]).Px()*NewAxis2_Px + (Object[i]).Py()*NewAxis2_Py +
                                          (Object[i]).Pz()*NewAxis2_Pz)/(Object[i]).P());
	 
	 
              mass1 *= NewAxis1_E/((NewAxis1_E+(Object[i]).E())*(NewAxis1_E+(Object[i]).E()));
	 
              mass2 *= NewAxis2_E/((NewAxis2_E+(Object[i]).E())*(NewAxis2_E+(Object[i]).E()));
	
	 
            if(mass1 < mass2) {
              if (Object_Group[i] != 1){ 
                I_Move = true;
              }
              Object_Group[i] = 1;
       
              Sum1_Px += (Object[i]).Px();
              Sum1_Py += (Object[i]).Py();
              Sum1_Pz += (Object[i]).Pz();
              Sum1_P += (Object[i]).P();
              Sum1_E += (Object[i]).E(); 
            } else {
              if (Object_Group[i] != 2){ 
                I_Move = true;
              }
              Object_Group[i] = 2;
              Sum2_Px += (Object[i]).Px();
              Sum2_Py += (Object[i]).Py();
              Sum2_Pz += (Object[i]).Pz();
              Sum2_P += (Object[i]).P();
              Sum2_E += (Object[i]).E(); 
	 
            }
      
      
          } else {
	
            if (Object_Group[i] == 1){
              Sum1_Px += (Object[i]).Px();
              Sum1_Py += (Object[i]).Py();
              Sum1_Pz += (Object[i]).Pz();
              Sum1_P += (Object[i]).P();
              Sum1_E += (Object[i]).E(); 
            }
            if (Object_Group[i] == 2){
              Sum2_Px += (Object[i]).Px();
              Sum2_Py += (Object[i]).Py();
              Sum2_Pz += (Object[i]).Pz();
              Sum2_P += (Object[i]).P();
              Sum2_E += (Object[i]).E(); 
            }
         
	
	
          }
	
	
        }
      }
    
    // recomputing the axes     

    Axis1[3] = sqrt(Sum1_Px*Sum1_Px + Sum1_Py*Sum1_Py + Sum1_Pz*Sum1_Pz);
    if (Axis1[3]<0.0001) {
        std::cout << "ZERO objects in group 1! " << std::endl; 
    } else {    
      Axis1[0] = Sum1_Px / Axis1[3];   
      Axis1[1] = Sum1_Py / Axis1[3];
      Axis1[2] = Sum1_Pz / Axis1[3];
      Axis1[4] = Sum1_E; 
    }
    
   
    
    Axis2[3] = sqrt(Sum2_Px*Sum2_Px + Sum2_Py*Sum2_Py + Sum2_Pz*Sum2_Pz);
    if (Axis2[3]<0.0001) {
        std::cout << "ZERO objects in group 2!";
    } else {  
      Axis2[0] = Sum2_Px / Axis2[3];   
      Axis2[1] = Sum2_Py / Axis2[3];
      Axis2[2] = Sum2_Pz / Axis2[3];
      Axis2[4] = Sum2_E; 
    }

    
  }
    return Object_Group;
}

package arrays;

import java.util.ArrayList;

public class AtMostKDuplicates {

	static int removeDuplicates(ArrayList<Integer> a, int at_most) {
	    //Here at most duplicates allowed is 2, can be scaled to N duplicates
	    //int at_most = 2;
	    //Index pointer to new modified array.
	    int j = 0;
	    //Handle corner cases of size less than 2.
	    if(a.size() <=2)
	        return a.size();
	        
	    for(int i=0; i<a.size(); i++)
	    {   
	        //Since corner cases are handled above, blindly add first two elements into same
	        //modified array.
	        //After filling first two elements check if the current element > 2 elements before the added array element.
	        if(j < at_most || a.get(i) > a.get(j-at_most))
	        {
	            a.set(j++,a.get(i));
	        }
	    }
	    System.out.println(a);
	    return j;
		
//		int prev = a.get(0);
//        int count = 1;
//        int j=1;
//       //  System.out.println(a.size());
//       //  System.out.println("after remove");
//       //  a.remove(0);
//       //  System.out.println(a);
//       //  System.out.println(a.size());
//       for(int i=1;i<a.size();i++)
//       {
//           if(a.get(i) == prev)
//           {
//               if(count < 2)
//               {
//                   count++;
//                   j++;
//               }
//               else
//               {
//                   //System.out.println(i + "after remove");
//                   a.remove(i);
//                   //System.out.println(a);
//                   i--;
//                   
//               }
//           }
//           else
//           {
//               prev = a.get(i);
//               count = 1;
//               j++;
//           }
//           
//           // if(j<2 || a.get(i) > a.get(j-2))
//           //      a.set(j++,a.get(i));
//       }
//       System.out.println(a);
//        return j; 
	}

	public static void main(String[] args) {
		ArrayList<Integer> a = new ArrayList<>();
		a.add(1);
		a.add(1);
		a.add(2);
		a.add(2);
		a.add(2);
		a.add(2);
		a.add(3);
		System.out.println(removeDuplicates(a, 2));
		

	}

}

package arrays;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

//Given a list of non negative integers, arrange them such that they form the largest number.
//
//For example:
//
//Given [3, 30, 34, 5, 9], the largest formed number is 9534330.
//
//Note: The result may be very large, so you need to return a string instead of an integer.
public class LargestNumber {
	
	 static String largestNumber(final List<Integer> a) {
	        if(a.size() == 0)
	            return "";
	         int n = a.size();  
	        String sa[] = new String[n];
	        
	        for(int i=0;i<n;i++)
	        {
	            sa[i] = Integer.toString(a.get(i));
	        }
	        
	        Arrays.sort(sa, new Comparator<String>(){
	           public int compare(String s1, String s2)
	           {
	               String val1 = s1+s2;
	               String val2 = s2+s1;
	               return val2.compareTo(val1);
	           }
	        });
	        
	        if(sa[0].equals("0"))
	        {
	            return "0";
	        }
	        
	        StringBuilder sb = new StringBuilder();
	        for(String s:sa)
	        {
	            sb.append(s);
	        }
	        return sb.toString();
	    }

	public static void main(String[] args) {
		//int a[] = {3, 30, 9, 5, 34};
		

	}

}

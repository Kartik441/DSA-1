package arrays;

import java.util.HashSet;
import java.util.Set;


public class ZeroSumSubarray {
	
	static boolean zeroSumSubarray(int []arr)
	{
		Set<Integer> set=new HashSet<Integer>();
		int sum=0;
		for(int element:arr)
		{
			set.add(sum);
			sum+=element;
			if(set.contains(sum))
				return true;
			
		}
		return false;
	}
	
	
	// can be solved using two pointer as well
	static boolean kSumSubarray(int []arr,int k)
	{
		Set<Integer> set=new HashSet<Integer>();
		int sum=0;
		for(int element:arr)
		{
			set.add(sum);
			sum+=element;
			if(set.contains(sum-k))
				return true;
			
		}
		return false;
	}

	public static void main(String[] args) {
		int arr[]= {2,2,3,-4,-2};
		int k=-3;
		System.out.println(zeroSumSubarray(arr));
		System.out.println(kSumSubarray(arr,k));

	}

}

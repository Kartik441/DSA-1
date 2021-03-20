package arrays;

public class LargestSumSubarray {
	
	static int kadanesAlgorithm(int a[])
	{
		// For an array consisting of all negative values, it will always give answer as 0 if handled incorrectly
		int n=a.length;
		int maxSum=Integer.MIN_VALUE;
		int curSum=0;
		for(int i=0;i<n;i++) {
			curSum+=a[i];
			
			if(curSum>maxSum)
				maxSum=curSum;
			
			if(curSum<0)
				curSum=0;
			
		}
		 System.out.println("Kadanes Sum: "+maxSum);
		 return maxSum;
	}
	
	
	// a variation could be to find largest sum subarray in circular array
	static void KadaneOnCircularArray(int a[])
	{
		// max subarray sum = Total sum of array - sum of non-contributing elements
		// we can use kadanes algorithm on same array with negative signs
		// to find sum of non contributing elements
		// it is giving 0 answer for all negative integers
		int n = a.length;
		int nonWrapSum = kadanesAlgorithm(a);
		int sum =0 ;
		for(int i=0;i<n;i++)
		{
			sum += a[i];
			//reversing the signs side wise
			a[i] = -a[i];
		}
		// now finding sum for non contributing elements
		
		int wrapSum = sum + kadanesAlgorithm(a);
		if(wrapSum == 0)  
			System.out.println(nonWrapSum);
		else
			System.out.println(Math.max(wrapSum, nonWrapSum));
		//it is giving answer as 0 for all negative integers
	}
	
	// Other Variations:
	// Maximum sum subarray after altering the array
	// Maximum sum subarray by flipping signs of at most k array elements 
	// Both present on gfg
	public static void main(String []args )
	{
		//KADEN'S ALGORITHM
		// For an array consisting of all negative values, it will always give answer as 0 if handled incorrectly
		int a[]= {-1, -3, -56, -34};//{-1,4,-2,-4,-1,3,5,-6,10};//{1,-3,2,5,7,-3,-2,10};
		kadanesAlgorithm(a);
		//int a1[] = {1, 4, -3, 8, -2, 5, 6};
		
		int a1[] = {4, -4, 6,-6, 10, -11, 12};
		int a2[] = {-4, 4, -6, 6, 10, 11, -12};
		kadanesAlgorithm(a1);
		kadanesAlgorithm(a2);
		KadaneOnCircularArray(a1);
		KadaneOnCircularArray(a);
//		    int []sum =new int[n];
//		    int maxSum=Integer.MIN_VALUE;
//		    int curSum=0;
//		    sum[0]=a[0];
//		    for(int i=1;i<n;i++)
//		    	sum[i]=sum[i-1]+a[i];
//		    for(int i=0;i<n;i++)
//		    {
//		    	for(int j=0;j<n;j++)
//		    	{
//		    		curSum=sum[j]-sum[i]+a[i];
//		    		if(curSum>maxSum)
//		    			maxSum=curSum;
//		    	}
//		    }
//		    System.out.println(maxSum);
	
		 
		 
		//BRUTEFORCE METHOD
		
//		int curSum=0;
//		int maxSum=Integer.MIN_VALUE;
//		for(int i=0;i<n;i++)
//		{
//			for(int j=0;j<n;j++)
//			{
//				curSum=0;
//				for(int k=i;k<=j;k++)
//				{
//					curSum+=a[k];
//				}
//				if(curSum>maxSum)
//					maxSum=curSum;
//					
//			}
//		}
//		System.out.println(maxSum);
	}


}

package arrays;

import java.util.Arrays;

// Given an array find the max j - i such that arr[j] > arr[i]
// One method is to use hashing to store the indexes and than sort the array
// Its time complexity: O(nlogn)


public class MaxIndexDiff {
	// another method is to traverse from the right of the array and keep track of max no to 
	// the right of current index including self
	// than we can use modified binary search

	// TC : O(nlogn)
	// SC : O(n)
	static void solve(int a[])
	{
		int n = a.length;
		int []maxFromEnd = new int[n];
		Arrays.fill(maxFromEnd, Integer.MIN_VALUE);
		
		maxFromEnd[n-1] = a[n-1]; 
		for(int i = n-2; i>=0 ; i--)
		{
			maxFromEnd[i] = Math.max(a[i], maxFromEnd[i+1]);
		}
		
		int res = -1;
		for(int i=0;i<n;i++)
		{
			// for every a[i], find the a[j] as per the conditions
			// a[j] can be find in O(logn) using binary search 
			int low = i+1;
			int high = n-1;
			int ans = -1;
			while(low <= high)
			{
				int mid = (low + high)/2;
				if(a[i] <= maxFromEnd[mid]) 
				{
					ans = mid - i;
					low = mid+1;
				}
				else
				{
					high = mid-1;
				}
			}
			
			res = Math.max(ans, res);
		}
		System.out.println(res);
		
		
	}
	
	
	
	
// TC: O(n) 
// SC : O(n)	
	static int maxIndexDiff(int a[])
	{
		int n = a.length;
		int i, j;
		
		int RMax[] = new int[n];
		int LMin[] = new int[n];
		
		LMin[0] = a[0];
		for(i = 1;i<n;i++)
		{
			LMin[i] = Math.min(a[i], LMin[i-1]);
		}
		RMax[n-1] = a[n-1];
		for(j=n-2;j>=0;j--)
		{
			RMax[j] = Math.max(a[j], RMax[j+1]);
		}
		
		// traverse both arrays from left to right to find optimal j-i 
		// Two pointer approach
		
		i = 0;
		j = 0;
		int maxDiff = -1;
		while(j<n && i<n)
		{
			if(LMin[i] < RMax[j])
			{
				// this can be the answer
				maxDiff = Math.max(maxDiff, j-i);
				// there is still chance that moving j might maximize the difference
				j = j+1;
			}
			else
			{
				// condition is not matching
				// no point in incrementing j bcoz the value is either going to decrease or be same
				// increment i bcoz the value will only increase or remain same which might get to the 
				// condition for the answer
				i = i+1;
			}
		}
		return maxDiff ;
	}
	
	public static void main(String[] args) {
	
		int a[] = {9, 2, 3, 4, 5, 6, 7, 8, 18, 0};
		int a1[] = {1,33, 5, 23, 4, 54, 55};
		int a2[] = {34, 8, 10, 3, 2, 80, 30, 33, 1};
		int a3[] = {6, 5, 4, 3, 2, 1};     
		
		System.out.println(maxIndexDiff(a));
		solve(a);
		System.out.println(maxIndexDiff(a1));
		solve(a1);
		System.out.println(maxIndexDiff(a2));
		solve(a2);
		System.out.println(maxIndexDiff(a3));
		solve(a3);

	}

}

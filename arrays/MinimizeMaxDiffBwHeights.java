package arrays;

import java.util.Arrays;
// better solution is available on stackoverflow by taking average of heights
public class MinimizeMaxDiffBwHeights {
	
	static int solve(int a[], int k)
	{
		int n = a.length;
		if(n == 1)
			return 0;
		Arrays.sort(a);
		
		// initializing result IMP: for two elements in array
		int ans = a[n-1] - a[0];
		
		int small = a[0] + k;
		int big = a[n-1] - k;
		
		// if after modifications values changed
		if(small > big)
		{
			int temp = small;
			small = big;
			big = temp;
		}
		
		for(int i=1;i<n-1;i++)
		{
			int subtract = a[i] - k;
			int add = a[i] + k;
			
			// if no change in small and big
			
			if(subtract >= small || add <= big)
				continue;
			
			if(subtract < small)
				small = subtract;
			else
				big = add;
		}
		return Math.min(ans, big - small);
		
	}

	public static void main(String[] args) {
		//int a[] = {1, 15, 10};
		//int a[] = {4, 6};
		//int a[] = {6, 10};
		//int a[] = {1, 10, 14, 14, 14, 15};
		int a[] = {1,2, 3};
		int k = 2;
		
		System.out.println(solve(a, k));

	}

}

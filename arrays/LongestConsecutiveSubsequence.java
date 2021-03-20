package arrays;

import java.util.HashMap;
import java.util.HashSet;
// This function returns the value of the size of longest consecutive subsequence found
public class LongestConsecutiveSubsequence {
	// this will not work for duplicates
	static int LCS(int []arr)
	{
		int n=arr.length;
		HashSet<Integer> s = new HashSet<Integer>(); 
		int finalans=0;
		
		for(int i=0;i<n;i++)
		{
			s.add(arr[i]);
		}
		for(int i=0;i<n;i++)
		{
			int l=0;
			if(!s.contains(arr[i]-1))
			{
				int e=arr[i];               //This might be the first element of subsequence
				while(s.contains(e))
				{
					e++;
					l++;
				}
			}
			if(l>finalans)
				finalans=l;
		}
		return finalans;
	}
	
	static int LCSD(int []a)
	{
		int n=a.length;
		HashMap<Integer, Integer>hm = new HashMap<>(); 
		int finalans=0;
		// using hashmap instead of hashset to handle duplicates
		for(int i=0;i<n;i++)
		{
			if(hm.containsKey(a[i]))
			{
				hm.put(a[i], hm.get(a[i])+1);
			}
			else
				hm.put(a[i], 1);
		}
		for(int i=0;i<n;i++)
		{
			int l=0;
			if(!hm.containsKey(a[i]-1))
			{
				int e=a[i];               //This might be the first element of subsequence
				while(hm.containsKey(e))
				{
					l += hm.get(e);
					e++;
				}
			}
			if(l>finalans)
				finalans=l;
		}
		return finalans;
	}

	public static void main(String[] args) {
		int arr1[]= {1,9,3,10,4,20,2,5};
		int arr[]= {36,41,56,35,44,33,34,92,43,32,42,37};
		
		System.out.println(LCS(arr1));
		int a[] = {10, 9, 4, 5, 4, 8, 6};
		System.out.println(LCS(a));
		
		int a1[] = {1, 2, 3, 2, 3, 7, 2, 1};
		System.out.println(LCSD(arr1));
		System.out.println(LCSD(a));
		System.out.println(LCSD(a1));
}
}

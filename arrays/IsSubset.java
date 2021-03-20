package arrays;

import java.util.HashSet;

//1.Naive method:o(mn)
//2. Sorting and using binary Search
//O(mlogm + nlogm)
//3. Sort both the arrays and then use merging concept   This can handle duplicates
//O(mlogm + nlogn)
//4. Use hashing(implemented below), can't be used for duplicates
//O(n) and space complexity O(n)


public class IsSubset {
	
	static boolean isSubset(int a1[], int a2[])
	{
		int n1 = a1.length;
		int n2 = a2.length;
		
		HashSet<Integer> hs = new HashSet<Integer>();
		
		for(int i=0; i<n1; i++)
		{
			if(!hs.contains(a1[i]))
				hs.add(a1[i]);
		}
		
		for(int i=0; i<n2; i++)
		{
			if(!hs.contains(a2[i]))
				return false;
		}
		
		return true;
	}

	public static void main(String[] args) {
		
		int arr1[] = {11, 1, 13, 21, 3, 7};
		int arr2[] = {11, 3, 7, 1};
		
		System.out.println(isSubset(arr1, arr2));
		

	}

}

package backtracking;

import java.util.ArrayList;

public class CombinationalSum {
	// given an array of positive integers arr[] and a sum x, find all unique combinations 
	// in arr where sum is equal to x
	// same as total ways coin change problem(combination), here we have to print combinations
	
	static ArrayList<ArrayList<Integer>> solve(int a[], int i, int sum, ArrayList<ArrayList<Integer>> ansList, ArrayList<Integer> ans)
	{
		if(i >= a.length)
			return ansList;
		
		if(sum < 0)
			return ansList;
		
		if(sum == 0)
		{
			// IMPORTANT
			// it is needed to make deep copy of ans as while backtracking the elements will be removed
			// and final list will contain empty arrayList
			// To make deep copy of array this method can be used
			//Arrays.copyOf(ans, ans.size())
			
			
			ArrayList<Integer> copy = new ArrayList<>();
			copy.addAll(ans);
			ansList.add(copy);			
			return ansList;
		}
		
		// one choice is to add the element
		ans.add(a[i]);
		solve(a, i, sum - a[i], ansList, ans);
		// backtrack
		ans.remove(ans.size()-1);
		
		// other choice is to not add the element
		solve(a, i+1, sum, ansList, ans);
		
		return ansList;
		
	}
	
	static void solve(int a[], int sum)
	{
		ArrayList<ArrayList<Integer>> ansList = new ArrayList<>();
		ArrayList<Integer> ans = new ArrayList<>();
		
		System.out.println(solve(a, 0, sum, ansList, ans));
	}

	public static void main(String[] args) {
		
		int a[] = {2, 4, 6, 8};
		solve(a,8);
	}

}

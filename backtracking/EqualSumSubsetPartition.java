package backtracking;

import java.util.ArrayList;

// partition an array into k subsets with equal sum

public class EqualSumSubsetPartition {
	
	static void solve(int a[], int k, int pt, int sets, ArrayList<ArrayList<Integer>> ans,int sum[] )
	{
		// some important base cases to cover all test cases
		// if k == 1 than whole array is your answer
		
		// if (k > n || sum % k != 0) than not possible
		
		
		if(pt == a.length)
		{		
			if(sets == k)
			{
				// check for equal sum and directly print than becomes answer of
				// K subset partition
				boolean flag = true;
				for(int i=1;i<sum.length;i++)
				{
					if(sum[i] != sum[i-1])
						flag = false;
				}
				if(flag)
					System.out.println("Answer is: "+ans);
			}
			return;
		}
		
		for(int i=0;i<k;i++)
		{
			
			if(ans.get(i).size() > 0)
			{
				
				ans.get(i).add(a[pt]);				
				sum[i] += a[pt];
				solve(a, k, pt+1, sets, ans, sum);
				// backtrack
				sum[i] -= a[pt];
				ans.get(i).remove(ans.get(i).size() - 1);
				
			}
			else
			{
				
				ans.get(i).add(a[pt]);
				
				sum[i] += a[pt];
				solve(a, k, pt+1, sets+1, ans, sum);
				// backtrack
				sum[i] -= a[pt];
				
				ans.get(i).remove(ans.get(i).size() - 1);
				
				break;          // this break statement is very important to avoid permutations
				
			}
			
		}
	}
	
	static void solve(int []a,int k)
	{
		ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
		for(int i=0;i<k;i++)
		{
			ans.add(new ArrayList<Integer>());
		}
		int sum[] = new int[k];
		solve(a, k, 0, 0, ans, sum);
	}

	public static void main(String[] args) {
		int a[] = {2, 1, 4, 5, 6};
		solve(a, 3);
		
		int a1[] = {2, 1, 5, 5, 6};
		solve(a1, 3);
		

	}

}

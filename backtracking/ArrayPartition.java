package backtracking;
// partition an array into two sets of equal sums
import java.util.ArrayList;

public class ArrayPartition {
	
	static boolean partition(int a[])
	{
		int sum = 0;
		for(int i=0;i<a.length;i++)
		{
			sum+= a[i];
		}
		//System.out.println(sum);
		ArrayList<Integer> ans = new ArrayList<Integer>();
		
	    if((sum&1) == 1)// for checking even
	    {
	    	System.out.println("Not possible");
	    	return false;
	    }
	    else
	    {
	    	if(partition(a, 0, sum/2, ans))
	    	{
	    		System.out.println(ans);
	    		return true;
	    	}
	    	else
	    		return false;
	    }
	}
	
	static boolean partition(int a[], int i, int sum, ArrayList<Integer> ans)
	{
		if(i >= a.length || sum<0)
			return false;
		
		if(sum == 0)
			return true;
		
		ans.add(a[i]);
		sum -= a[i];
		
		boolean isPossible = partition(a, i+1, sum, ans);
		
		if(isPossible)
			return true;
		
		ans.remove(ans.size() - 1);
		sum += a[i];
		
		return partition(a, i+1, sum, ans);
		
	}

	public static void main(String[] args) {
		int a[] = {2, 1, 2, 3, 4, 8};
		int a1[] = {2, 1, 5, 2};
		
		System.out.println(partition(a));
		System.out.println(partition(a1));
		

	}

}

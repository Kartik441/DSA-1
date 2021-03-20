package arrays;
// Maximum difference between two elements such that greater element comes after smaller one

public class MaxDiff {
// This can solve stock problem if we are allowed to buy and sell only once	
	static int maxDiff(int a[])
	{
		int n = a.length;
		int max_diff = a[1] - a[0];
		int min_element = a[0];
		for(int i=1;i<n;i++)
		{
			if(a[i] - min_element  > max_diff)
				max_diff = a[i] - min_element;
			
			if(a[i] < min_element)
				min_element = a[i];
		}
		return max_diff;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int a[] = {1, 2, 90, 10, 110};
		int a1[] = {2, 3, 10, 6, 4, 8, 1};
		int a2[] = {7, 9, 5, 6, 3, 2};
		int a3[] = {6, 5, 4, 3, 2, 1};
		int a4[] = {9, 2, 3, 4, 5, 6, 7, 8, 18, 0};
		System.out.println((maxDiff(a4)));

	}

}

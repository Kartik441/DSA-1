package arrays;

//1 Method: Naive method will take O(mn)
//2 Method: Use binary search  O(mlogn)
//3 Method: Optimized bs (implemented below)



public class Max1s {
	
	static int findRow(int a[][])
	{
		int m = a[0].length;        // no. of columns
		int n = a.length;          // no. of rows
		int max = Integer.MIN_VALUE;
		
		int row = -1;
		
		for(int i=0; i<n; i++)
		{
			int ones = bsAsset(a[i]);
			
			
			if(ones > max) {
				max = ones;
				row = i;
			}
		}
		
		return row;
		
	}
	
	static int bsAsset(int a[])
	{
		int i = 0;
		int j = a.length - 1;
		int ans = a.length;
		
		while(i<=j)     // WARNING: not (i<j)
		{
			int mid = (i+j)/2;
		
			if( a[mid] == 1) {
				ans = mid;
				j = mid -1;
			}
			
			else
			{
				i = mid+1;
			}
			
			
		}
		
		
		
		return a.length - ans;
		
		
	}

	public static void main(String[] args) {
		int a[][] = {
				{0, 1, 1, 1},
				{0, 0, 1, 1},
				{0, 1, 1, 1},
				{0, 0, 0, 1},
				{1, 1, 1, 1},
				{0, 0, 0, 0}
		
		};
		
		System.out.println(findRow(a));
		
		
	}

}

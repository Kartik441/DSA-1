package arrays;

public class ThreeWayPartition {
// segregate 0,1, and 2
// this is dutch flag algorithm
	static void segregate(int a[])
	{
		int n = a.length;
		int l=0;
		int m=0;
		int h=n-1;
		
		while(m<=h)
		{
			if(a[m] == 0)
			{
				// swap a[l] and a[m] and increment l and m
				int temp = a[l];
				a[l] = a[m];
				a[m] = temp;
				
				l++;
				m++;
			}
			
			if(m>=h)
				break;
			if(a[m] == 1)
			{
				m++;
			}
			if(m>=h)
				break;
			if(a[m] == 2)
			{
				// swap a[m] and a[h] and decrement h
				int temp = a[m];
				a[m] = a[h];
				a[h] = temp;
				
				h--;
			}
		}
		for(int i=0;i<n;i++)
		{
			System.out.print(a[i] + " ");
		}
		System.out.println();
		
	}
	
	// we can use the above algorithm for solving 3 way partitioning 
	static void solve(int a[], int low, int high)
	{
		int n = a.length;
		int l=0;
		int m=0;
		int h=n-1;
		
		while(m<=h)
		{
			if(a[m] < low)
			{
				// swap a[l] and a[m]
				int temp = a[l];
				a[l] = a[m];
				a[m] = temp;
				l++;
				m++;
			}
			if(m >= h)
				break;
			if(a[m] >= low && a[m] <= high)
			{
				m++;
			}
			if(m >= h)
				break;
			if(a[m] > high)
			{
				// swap a[m] and a[h]
				int temp = a[m];
				a[m] = a[h];
				a[h] = temp;
				h--;
			}
		}
		for(int i=0;i<n;i++)
		{
			System.out.print(a[i] + " ");
		}
		System.out.println();
	}
	public static void main(String[] args) {
	
		int a[] = {2, 2, 1, 1, 0, 2, 2, 1, 0,2};
		segregate(a);
		
		int a1[] = {1, 14, 5, 20, 4, 2, 54, 20, 87, 98, 3, 1, 32};
		solve(a1, 14, 20);
		
		int a2[] = {1, 14, 5, 20, 4, 2, 54,20, 87, 98, 3,1, 32};
		solve(a2, 20, 20);

	}

}

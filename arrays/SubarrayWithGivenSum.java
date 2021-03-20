package arrays;
// this can be done using hashing as well (it will handle negative numbers as well)
// here is using two pointer algorithm (sliding window) 
public class SubarrayWithGivenSum {
// also known as k sum subarray
	
	// this function will give only first subarray with sum equal to k
	static void solve(int a[], int sum)
	{
		int n = a.length;
		
		int i=0, j=0;
		int st = -1, end = -1, currSum = 0;
		
		while(j<n && currSum + a[j] <= sum)
		{
			currSum += a[j];
			j++;
		}
		
		if(currSum == sum)
		{
			System.out.println(i+1 + " "+j);
			return;
		}
		while(j<n)
		{
			currSum += a[j];
			while(currSum > sum)
			{
				currSum -= a[i];
				i++;
			}
			if(currSum == sum)
			{
				st = i+1;
				end = j+1;
				break;
			}
			j++;
		}
		System.out.println(st+" "+end);
		
	}
	// this function will give all subarrays with sum equal to k
	static void solve1(int a[], int k)
	{
		int n = a.length;
		int i = 0;
		int j = 0;
		int sum = a[0];
		while(i < n && j<n)
		{
			// expanding the window
			while(sum < k)
			{
				j++;
				if(j<n)
				{
					sum += a[j];
				}
				else 
					break;
			}
			
			if(sum == k)
			{
				System.out.println("start & end " + i+" "+j);
				// new starting point for another answer if possible
				sum -= a[i];
				i++;
			}
			
			while(sum > k && i <= j)
			{
				sum -= a[i];
				i++;
			}
			
			if(sum == k)
			{
				System.out.println("start & end " + i+" "+j);
				sum -= a[i];
				i++;
			}
			
		}
		
	}

	public static void main(String[] args) {
		
		int a[] = {1, 2, 3, 7, 5};
		solve(a, 12);
		solve1(a, 12);
		int a1[] = {1, 2, 3, 8};
		solve(a1, 5);
		solve1(a1, 5);
		
		int a2[] = {1, 2, 34, 5};
		solve(a2, 34);
		solve1(a2, 34);

	}

}

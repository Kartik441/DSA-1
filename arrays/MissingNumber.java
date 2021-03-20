package arrays;

public class MissingNumber {
	
	//list of n-1 integers ranging form 1 to n.No duplicates.Find missing no.
	static int missingNumber(int arr[])
	{
		int sum=0;
		int n=arr.length+1;
		int tsum=(n*(n+1))/2;
		for(int i:arr)
			sum+=i;
		return tsum-sum;
		
	}
	//Another Method:Using XOR in bit manipulation package, which is faster 
	// idea is to create duplicates, than all the numbers which are already present in the
	// array will cancel out each other and we will be left with answer
	static void missingXor(int a[])
	{
		int n = a.length;
		int xor = 0;
		for(int i=0;i<n;i++)
		{
			xor = xor^a[i];
		}
		
		// creating duplicates
		for(int i=1;i<=n+1;i++)
		{
			xor = xor^i;
		}
		System.out.println(xor);
	}

	public static void main(String[] args) {
		int arr[]= {1,2,4,5,6,3,8,9};
		System.out.println(missingNumber(arr));
		missingXor(arr);

	}

}

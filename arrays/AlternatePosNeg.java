package arrays;
// Variation of this could be to maintain order while arranging
public class AlternatePosNeg {
	
	static int[] alternateUnordered(int []arr)
	{
		int i=0;
		int j=arr.length-1;
		
		while(i<=j)
		{
			while(i<arr.length &&arr[i]<0)
				i++;
			while(j>=0 && arr[j]>0)
				j--;
			
			// swapping values at i and j
			if(i<=j) {      // This condition is important
				int temp=arr[i];
				arr[i]=arr[j];
				arr[j]=temp;
			}
			
		}
		return arr;
	}
	
	// 2 part
	// maintain order while rearranging
	
	// 1 Method : Modified insertion sort
	// TC: O(n^2)
	// SC : O(1)
	
	// 2 Method
	// Ln and Lp refers to negative and positive part of left subarray
	// Rn and Rp refers to negative and positive part of right subarray
	
	// 1. reverse Lp and Rn separately. 
	// [Ln Lp Rn Rp] -> [Ln Lp'Rn'Rp]
	// 2. reverse [Lp'Rn'] 
	// [Ln Lp'Rn'Rp] -> [Ln Rn Lp Rp]
	
	static void rearrange(int a[], int l, int r)
	{
		if(l < r)
		{
			int mid = (l+r)/2;
			
			rearrange(a, l, mid);
			rearrange(a, mid+1, r);
			
			merge(a, l, mid, r);
		}
	}
	
	static void merge(int a[], int l, int m, int r)
	{
		int i = l;
		int j = m + 1;
		
		while(i<=m && a[i] < 0)    //move until hit first positive number
			i++;
		
		// a[i...m] is positive
		
		while(j<=r && a[j] < 0)
			j++;
		
		// a[j..r] is positive
		// now apply the merge approach
		// 1 step
		reverse(a, i, m);
		reverse(a, m+1, j-1);
		
		// 2 step
		reverse(a, i, j-1);
	}
	
	static void reverse(int arr[],int a,int b)
	{
	
		int i=a;                                                      
		int j=b;
		while(i<j)
		{
			arr[i]=arr[i]+arr[j];
			arr[j]=arr[i]-arr[j];
			arr[i]=arr[i]-arr[j];
			i++;
			j--;
		}
		
		
	}

	public static void main(String[] args) {
		int arr[]= {-1,3,5,-6,245,-56,78,45,-23,-58,27,-7,8};
		arr=alternateUnordered(arr);
		for(int i:arr)
			System.out.print(i+" ");
		System.out.println();
		int arr1[]= {-1,3,5,-6,245,-56,78,45,-23,-58,27,-7,8};
		rearrange(arr1, 0, arr1.length-1);
		for(int i:arr1)
			System.out.print(i+" ");


	}

}

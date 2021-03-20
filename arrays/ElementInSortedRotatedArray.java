package arrays;

public class ElementInSortedRotatedArray {
	
	static int minimumElement(int []arr)
	{

		int n=arr.length;
		int low=0;
		int high=n-1;
		while(low<high)
		{
			int mid=(low+high)/2;
			if(arr[mid]==arr[high])
				high--;
			else if(arr[mid]>arr[high])
				low=mid+1;
			else
				high=mid;       // NOT mid-1
		}
		return arr[high];
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int arr[]= {5,6,1,2,3,4};
		System.out.println(minimumElement(arr));

	}

}

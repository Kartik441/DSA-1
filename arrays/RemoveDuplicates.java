package arrays;

public class RemoveDuplicates {
	// In sorted order
	static int[] removeDuplicates(int arr[])
	{
		int i=0;
		int j=0;
		for(i=0;i<arr.length-1;i++)
		{
			if(arr[i]!=arr[i+1])
			{
				arr[j]=arr[i];
				j++;
			}
		}
		arr[j]=arr[i];
		int arr1[]=new int[j+1];
		for(int k=0;k<=j;k++)
		{
			arr1[k]=arr[k];
			
		}
		arr=null;
		
		return arr1;
	}

	public static void main(String[] args) {
		
		int arr[]= {1,2,2,3,3,3,3,4,5,5,5,5,5,5,5,6};
		arr=removeDuplicates(arr);
		for(int i:arr)
			System.out.print(i+" ");
	}
}
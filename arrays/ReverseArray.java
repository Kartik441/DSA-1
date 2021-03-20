package arrays;

public class ReverseArray {
	
	static int[] reverseArray(int arr[],int a,int b)
	{
		int n=arr.length;
		int i=a;                                                           //  0;
		int j=b;
		while(i<j)
		{
			arr[i]=arr[i]+arr[j];
			arr[j]=arr[i]-arr[j];
			arr[i]=arr[i]-arr[j];
			i++;
			j--;
		}
		
		return arr;
	}

	public static void main(String[] args) {
		int arr[]= {1,2,3,4,5,6,7,8,9};
		
		for(int i:arr)
			System.out.print(i+" ");
		System.out.println();
		arr=reverseArray(arr,2,6);                         
		for(int i:arr)
			System.out.print(i+" ");
		System.out.println();

	}

}

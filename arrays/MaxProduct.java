package arrays;

public class MaxProduct {
	
	static int maxProduct(int []arr)
	{
		int n=arr.length;
		int maxVal=arr[0];
		int minVal=arr[0];
		
		int maxProduct = arr[0];
		
		for(int i=1;i<n;i++)
		{
			if(arr[i]<0)
			{
				int temp=maxVal;
				maxVal=minVal;
				minVal=temp;
			}
			
			maxVal = Math.max(arr[i], maxVal*arr[i]);
			minVal = Math.min(arr[i],minVal*arr[i]);
			
			System.out.println("MIN: "+ minVal);
			System.out.println("MAX: "+ maxVal);
			System.out.println("MAXPRO: "+ maxProduct);
			
			maxProduct = Math.max(maxProduct, maxVal);
			System.out.println("MAXPRO: "+ maxProduct);
		}
		
		return maxProduct;
		
	}

	public static void main(String[] args) {
		
		int arr[]= {-150,-3,-1,-5};       
		int a1[] = {-4, 0, -5, 0};
		System.out.println(maxProduct(a1));
		

	}

}

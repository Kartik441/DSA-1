//return true if the sum of three elements is equal to zero   

package arrays;

import java.util.Arrays;

public class TripleSum {
	public static boolean twoSum(int []a,int x,int i)
	{
		int j=a.length-1;
		while(i<j)
		{
			if(a[i]+a[j]<x)
				i++;
			else if(a[i]+a[j]>x)
				j--;
			else 
				return true;
		}
		return false;
	}

	public static void main(String[] args) {
		int[] a= {1,-3,-2,4,6,8,11,-10};
		int flag=0;
		Arrays.sort(a);
		for(int i=0;i<a.length;i++)
		{
			if(twoSum(a,-a[i],i+1))
			{
				System.out.println("True");
				flag=1;
				break;
			}
			
		}
		if(flag==0)
		{
			System.out.println("false");
		}
		

	}

}

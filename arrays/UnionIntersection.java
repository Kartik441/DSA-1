package arrays;
import java.util.*;

public class UnionIntersection {
	public static Vector<Integer> union(int []a,int []b)
	{
		Vector<Integer> union=new Vector<Integer>();
		int n=a.length;
		int m=b.length;
		int i=0,j=0;
		while(i<n&&j<m)
		{
			if(a[i]<b[j])
			{
				union.add(a[i]);
				i++;
			}
			else if(a[i]>b[j])
			{
				union.add(b[j]);
				j++;
			}
			else   // Means both are equal
			{
				union.add(a[i]);
				i++;
				j++;
			}
	
		}
	
		while(i<n)
		{
			
				union.add(a[i]);
				i++;
		}
		while(j<m)
		{
				union.add(b[j]);
				j++;
			
		}
		return union;
	}
	public static Vector<Integer> intersection(int []a,int []b)
	{
		Vector<Integer> intersection =new Vector<Integer>();
		int n=a.length;
		int m=b.length;
		int i=0,j=0;
		while(i<n&&j<m)
		{
			if(a[i]<b[j])
			{
		        i++;
			}
			else if(a[i]>b[j])
			{
				j++;
			}
			else
			{
				intersection.add(a[i]);
				i++;
				j++;
			}
		}
		return intersection;
	}

public static  void main(String[] args)
{
	int a[]= {100,150,175,252,300,500};
	int b[]= {100,26,38,47,60,70,252,500};
	Vector<Integer> union=new Vector<Integer>(union(a,b));
	Vector<Integer> inter=new Vector<Integer>(intersection(a,b));
	System.out.println(union);
	System.out.println(inter);
}
	
}

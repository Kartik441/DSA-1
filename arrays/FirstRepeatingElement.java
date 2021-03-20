package arrays;

//n is the total no of elements ranging between 1 to n
public class FirstRepeatingElement {
// if we are not allowed to change array than we can use other index array	
	static int firstRepeatingElement(int arr[])
	{
		for(int i=0;i<arr.length;i++)
		{
			if(arr[Math.abs(arr[i])-1]>=0)      // if visiting first time
				arr[Math.abs(arr[i])-1]*=-1;    // add marker(make it negative)
			else
				return arr[i];        
		}
		return -1;
	}

	public static void main(String[] args) {
		int arr[]= {1,2,3,4,5,7,6,5,4,7};
		System.out.println(firstRepeatingElement(arr));

	}

}

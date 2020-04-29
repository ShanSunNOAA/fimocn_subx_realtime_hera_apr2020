module stringlist

  implicit none

  private

  public ::              &
    stringlist_add,      &
    stringlist_at_index, &
    stringlist_destroy,  &
    stringlist_find,     &
    stringlist_first,    &
    stringlist_has,      &
    stringlist_index_of, &
    stringlist_last,     &
    stringlist_node,     &
    stringlist_rest

  type stringlist_node
    character,pointer::chararr(:)
    type(stringlist_node),pointer::next=>null()
  end type stringlist_node

contains

  ! public

  subroutine stringlist_add(list,string)
    type(stringlist_node),pointer::list
    character(len=*),intent(in)::string
    type(stringlist_node),pointer::ptr
    ptr=>stringlist_blank(list)
    call stringlist_set(ptr,string)
  end subroutine stringlist_add

  function stringlist_at_index(list,index)
    type(stringlist_node),pointer::stringlist_at_index,list
    integer,intent(in)::index
    integer::i
    stringlist_at_index=>null()
    if (.not.associated(list)) return
    stringlist_at_index=>list
    i=1
    do while (associated(stringlist_at_index))
      if (i.eq.index) return
      stringlist_at_index=>stringlist_at_index%next
      i=i+1
    end do
    stringlist_at_index=>null()
  end function stringlist_at_index

  subroutine stringlist_destroy(list)
    type(stringlist_node),pointer::list,ptr
    do while (associated(list))
      ptr=>list%next
      deallocate(list)
      list=>ptr
    end do
  end subroutine stringlist_destroy

  function stringlist_find(list,string)
    type(stringlist_node),pointer::stringlist_find,list
    character(len=*),intent(in)::string
    stringlist_find=>null()
    if (.not.associated(list)) return
    stringlist_find=>list
    do while (associated(stringlist_find))
      if (stringlist_match(stringlist_find,string)) return
      stringlist_find=>stringlist_find%next
    end do
  end function stringlist_find

  function stringlist_first(list)
    type(stringlist_node),pointer::stringlist_first,list
    stringlist_first=>null()
    if (.not.associated(list)) return
    stringlist_first=>list
  end function stringlist_first

  function stringlist_has(list,string)
    logical::stringlist_has
    type(stringlist_node),pointer::list
    character(len=*),intent(in)::string
    stringlist_has=.false.
    if (associated(stringlist_find(list,string))) stringlist_has=.true.
  end function stringlist_has

  function stringlist_index_of(list,string)
    type(stringlist_node),pointer::list,ptr
    integer::stringlist_index_of
    character(len=*)::string
    stringlist_index_of=0
    if (.not.associated(list)) return
    ptr=>list
    do while (associated(ptr))
      stringlist_index_of=stringlist_index_of+1
      if (stringlist_match(ptr,string)) return
      ptr=>ptr%next
    end do
    stringlist_index_of=0
  end function stringlist_index_of

  function stringlist_last(list)
    type(stringlist_node),pointer::stringlist_last,list
    stringlist_last=>null()
    if (.not.associated(list)) return
    stringlist_last=>list
    do while (associated(stringlist_last%next))
      stringlist_last=>stringlist_last%next
    end do
  end function stringlist_last

  function stringlist_rest(list)
    type(stringlist_node),pointer::stringlist_rest,list
    stringlist_rest=>null()
    if (.not.associated(list)) return
    stringlist_rest=>list%next
  end function stringlist_rest

  ! private
  
  function stringlist_blank(list)
    type(stringlist_node),pointer::stringlist_blank,list
    if (associated(list)) then
      stringlist_blank=>list
      do while (associated(stringlist_blank%next))
        stringlist_blank=>stringlist_blank%next
      end do
      allocate(stringlist_blank%next)
      stringlist_blank=>stringlist_blank%next
    else
      allocate(list)
      stringlist_blank=>list
    end if
  end function stringlist_blank

  function stringlist_match(ptr,string)
    logical::stringlist_match
    type(stringlist_node),pointer::ptr
    character(len=*)::string
    integer::i
    stringlist_match=.true.
    do i=1,len(string)
      if (ptr%chararr(i).ne.string(i:i)) then
        stringlist_match=.false.
        return
      end if
    end do
  end function stringlist_match

  subroutine stringlist_set(ptr,string)
    type(stringlist_node),pointer::ptr
    character(len=*)::string
    integer::i
    allocate(ptr%chararr(len(string)))
    do i=1,len(string)
      ptr%chararr(i)=string(i:i)
    end do
  end subroutine stringlist_set

end module stringlist
